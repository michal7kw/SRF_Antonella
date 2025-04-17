import os
import git
import hashlib
import datetime
import json
import shutil
import glob
from pathlib import Path
import tempfile
import fitz  # PyMuPDF
# Removed PIL import as it's not needed for PDFs
import sys

def extract_git_pdfs(repo_path, output_dir):
    """
    Extract PDFs from all branches and commits in a local Git repository.

    Args:
        repo_path: Path to the local Git repository
        output_dir: Directory to store the extracted PDFs
    """
    print(f"Starting PDF extraction from local repository at {repo_path}")

    # Validate repo path
    repo_path_obj = Path(repo_path)
    if not (repo_path_obj / '.git').is_dir():
        print(f"Error: '{repo_path}' is not a valid Git repository.")
        return None

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    try:
        repo = git.Repo(repo_path)
    except git.InvalidGitRepositoryError:
        print(f"Error: Invalid Git repository at {repo_path}")
        return None
    except Exception as e:
        print(f"Error opening repository: {e}")
        return None

    # Create results directory structure (named 'pdfs')
    results_dir = output_path / 'pdfs'
    results_dir.mkdir(exist_ok=True)

    # Create metadata directory
    metadata_dir = output_path / 'metadata'
    metadata_dir.mkdir(exist_ok=True)

    # Get all local branches
    branches = [head.name for head in repo.heads]
    # Also include remote branches if needed, requires fetching first
    try:
        repo.git.fetch('--all')
        remote_branches = []
        for ref in repo.remote().refs:
            if 'HEAD' not in ref.name:
                 branch_name_parts = ref.name.split('/')
                 if len(branch_name_parts) > 1:
                     local_branch_name = '/'.join(branch_name_parts[1:])
                     if local_branch_name not in branches:
                         remote_branches.append(ref.name)
                         branches.append(local_branch_name)
                 elif ref.name not in branches:
                     branches.append(ref.name)

        all_refs_to_check = [head.name for head in repo.heads] + remote_branches
        all_refs_to_check = list(set(all_refs_to_check))

    except git.GitCommandError as e:
        print(f"Warning: Could not fetch remote branches: {e}. Processing local branches only.")
        all_refs_to_check = [head.name for head in repo.heads]


    print(f"Found {len(all_refs_to_check)} refs (local branches and remote refs) to check.")

    # Dictionary to store unique PDFs
    unique_pdfs = {}

    # Track overall statistics
    stats = {
        'total_refs': len(all_refs_to_check),
        'processed_refs': 0,
        'total_commits': 0,
        'total_pdfs': 0, # Renamed from total_images
        'unique_pdfs': 0 # Renamed from unique_images
    }

    # Store current branch to return later
    try:
        original_branch = repo.active_branch.name
    except TypeError: # Detached HEAD state
        original_branch = repo.head.commit.hexsha
        print("Warning: Repository is in a detached HEAD state.")


    # Process each ref (branch or remote ref)
    processed_commits = set()

    for ref_name in all_refs_to_check:
        local_branch_name = ref_name.split('/')[-1]
        print(f"\nProcessing ref: {ref_name} (as {local_branch_name})")
        stats['processed_refs'] += 1

        safe_branch_name = local_branch_name.replace('/', '_').replace('\\', '_')
        branch_dir = results_dir / safe_branch_name
        branch_dir.mkdir(exist_ok=True)

        try:
            commits = list(repo.iter_commits(ref_name))
            print(f"  Found {len(commits)} commits for this ref")
        except git.GitCommandError as e:
            print(f"  Error getting commit history for {ref_name}: {e}")
            continue
        except Exception as e:
            print(f"  Unexpected error getting commits for {ref_name}: {e}")
            continue

        for commit in commits:
            commit_hash = commit.hexsha
            if commit_hash in processed_commits:
                continue

            processed_commits.add(commit_hash)
            stats['total_commits'] += 1

            short_hash = commit_hash[:8]
            commit_date = datetime.datetime.fromtimestamp(commit.committed_date)
            commit_date_str = commit_date.strftime('%Y%m%d_%H%M%S')

            print(f"  Processing commit: {short_hash} from {commit_date_str}")

            commit_dir = branch_dir / f"{commit_date_str}_{short_hash}"

            try:
                tree = commit.tree
            except Exception as e:
                print(f"    Error accessing tree for commit {short_hash}: {e}")
                continue

            pdf_extension = '.pdf' # Only look for PDFs
            commit_pdfs_found = 0

            for blob in tree.traverse():
                if blob.type == 'blob' and Path(blob.path).suffix.lower() == pdf_extension:
                    rel_path = Path(blob.path)
                    commit_pdfs_found += 1
                    stats['total_pdfs'] += 1

                    commit_dir.mkdir(exist_ok=True)

                    try:
                        content = blob.data_stream.read()
                        file_hash = hashlib.md5(content).hexdigest()
                    except Exception as e:
                        print(f"    Error reading blob data for {rel_path} in commit {short_hash}: {e}")
                        continue

                    dest_file = commit_dir / f"{file_hash}{rel_path.suffix}"
                    thumbnail_filename = f"{file_hash}.png"
                    thumbnail_dest_file = commit_dir / thumbnail_filename
                    thumbnail_rel_path_str = None # Default

                    # Determine relative path for thumbnail (used in metadata and potentially HTML)
                    # Need safe_branch_name here
                    thumbnail_rel_path = Path('pdfs') / safe_branch_name / commit_dir.name / thumbnail_filename
                    potential_thumbnail_rel_path_str = thumbnail_rel_path.as_posix()

                    # --- Write PDF and Generate Thumbnail (if needed) ---
                    pdf_written_in_this_iteration = False
                    if not dest_file.exists():
                        try:
                            with open(dest_file, 'wb') as f_out:
                                f_out.write(content)
                            pdf_written_in_this_iteration = True
                        except Exception as e:
                            print(f"    Error writing PDF file {dest_file}: {e}")
                            # If PDF write fails, we can't generate a thumbnail, so skip to next blob
                            continue

                    # Try generating thumbnail only if PDF exists (written now or previously) and thumbnail doesn't
                    if dest_file.exists() and not thumbnail_dest_file.exists():
                        try:
                            doc = fitz.open(dest_file)
                            if len(doc) > 0:
                                page = doc.load_page(0)
                                pix = page.get_pixmap(dpi=72) # Render at 72 DPI for smaller size
                                pix.save(thumbnail_dest_file)
                                thumbnail_rel_path_str = potential_thumbnail_rel_path_str # Store the relative path
                            doc.close()
                        except Exception as thumb_err:
                            print(f"      Warning: Could not generate thumbnail for {rel_path} ({dest_file}): {thumb_err}")
                            # Don't delete the PDF if thumbnail fails, just skip thumbnail

                    # --- Store Metadata ---
                    if file_hash not in unique_pdfs:
                        # Check if thumbnail exists (could be from a previous run)
                        current_thumbnail_path = potential_thumbnail_rel_path_str if thumbnail_dest_file.exists() else None
                        unique_pdfs[file_hash] = {
                            'hash': file_hash,
                            'original_path': str(rel_path),
                            'size_bytes': len(content),
                            'thumbnail_path': current_thumbnail_path, # Store path if it exists
                            'occurrences': []
                        }
                        stats['unique_pdfs'] += 1
                    # If PDF entry exists, ensure thumbnail path is updated if generated now
                    elif unique_pdfs[file_hash].get('thumbnail_path') is None and thumbnail_dest_file.exists():
                         unique_pdfs[file_hash]['thumbnail_path'] = potential_thumbnail_rel_path_str


                    # --- Add Occurrence Info ---
                    occurrence_info = {
                        'branch': local_branch_name,
                        'ref': ref_name,
                        'commit': commit_hash,
                        'commit_short': short_hash,
                        'date': commit_date.isoformat(),
                        'path': str(rel_path)
                    }

                    existing_occurrences = unique_pdfs[file_hash]['occurrences']
                    is_duplicate_occurrence = any(
                        occ['commit'] == commit_hash and occ['path'] == str(rel_path) and occ['ref'] == ref_name
                        for occ in existing_occurrences
                    )

                    if not is_duplicate_occurrence:
                         unique_pdfs[file_hash]['occurrences'].append(occurrence_info)

                    try:
                        with open(dest_file, 'wb') as f_out:
                            f_out.write(content)
                    except Exception as e:
                        print(f"    Error writing PDF file {dest_file}: {e}")

            if commit_pdfs_found > 0:
                 print(f"    Found {commit_pdfs_found} PDFs in commit {short_hash}")

    print("\nRestoring original repository state...")
    try:
        repo.git.checkout(original_branch)
        print(f"Checked out original ref: {original_branch}")
    except git.GitCommandError as e:
        print(f"Error checking out original ref '{original_branch}': {e}")
        print("Repository might be in an unclean state.")
    except Exception as e:
         print(f"Unexpected error restoring state: {e}")

    metadata_file = metadata_dir / 'pdfs.json' # Changed filename
    print(f"\nSaving metadata to {metadata_file}...")
    try:
        with open(metadata_file, 'w') as f:
            json.dump(unique_pdfs, f, indent=2)
    except Exception as e:
        print(f"Error saving metadata: {e}")

    html_index_file = output_path / 'index.html'
    print(f"Creating HTML index at {html_index_file}...")
    try:
        create_pdf_html_index(output_path, unique_pdfs, stats) # Changed function call
    except Exception as e:
        print(f"Error creating HTML index: {e}")

    print(f"\nExtraction complete!")
    print(f"Processed {stats['processed_refs']} refs and {stats['total_commits']} unique commits")
    print(f"Found {stats['unique_pdfs']} unique PDFs out of {stats['total_pdfs']} total PDF instances found across history") # Updated labels
    print(f"Results stored in: {output_path.resolve()}")

    return output_dir

# Renamed function and adapted for PDFs
def create_pdf_html_index(output_path, unique_pdfs, stats):
    """Create an HTML index page for PDFs with search functionality"""
    results_dir = output_path / 'pdfs' # Changed directory name
    html_file_path = output_path / 'index.html'

    with open(html_file_path, 'w', encoding='utf-8') as f:
        f.write('''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Git Repository PDF Search</title> <!-- Changed title -->
    <style>
        body { font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Oxygen, Ubuntu, Cantarell, "Open Sans", "Helvetica Neue", sans-serif; margin: 0; padding: 20px; background-color: #f8f9fa; color: #212529; }
        .container { max-width: 1600px; margin: 0 auto; }
        h1, h2 { color: #343a40; }
        .stats { background-color: #e9ecef; padding: 15px; border-radius: 8px; margin-bottom: 20px; border: 1px solid #dee2e6; }
        .controls { background-color: #fff; padding: 20px; border-radius: 8px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        .search, .filters { margin-bottom: 15px; }
        label { display: block; margin-bottom: 5px; font-weight: 500; color: #495057; }
        input[type="text"], select { width: 100%; padding: 10px; box-sizing: border-box; font-size: 14px; border: 1px solid #ced4da; border-radius: 4px; }
        input[type="text"]:focus, select:focus { border-color: #80bdff; outline: 0; box-shadow: 0 0 0 0.2rem rgba(0,123,255,.25); }
        .filter-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; }
        /* Changed grid class name */
        .pdf-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(250px, 1fr)); gap: 20px; }
        /* Changed item class name */
        .pdf-item { background-color: #fff; border: 1px solid #dee2e6; border-radius: 8px; overflow: hidden; box-shadow: 0 1px 3px rgba(0,0,0,0.05); transition: box-shadow 0.2s ease-in-out; display: flex; flex-direction: column; }
        .pdf-item:hover { box-shadow: 0 4px 8px rgba(0,0,0,0.1); }
        /* Added img style for thumbnail */
        .pdf-item img { display: block; margin: 0 auto; background: #f9f9f9; border-bottom: 1px solid #eee; }
        .pdf-link-area { padding: 10px; text-align: center; background-color: #e9ecef; border-top: 1px solid #dee2e6; min-height: auto; } /* Adjusted padding */
        .pdf-link-area a { font-weight: bold; word-break: break-all; font-size: 13px; }
        .pdf-info { padding: 15px; font-size: 13px; line-height: 1.5; flex-grow: 1; } /* Added flex-grow */
        .pdf-info p { margin: 0 0 8px 0; word-wrap: break-word; }
        .pdf-info strong { color: #495057; }
        .pdf-info .path { font-family: monospace; font-size: 12px; color: #6c757d; }
        .pdf-info a { color: #007bff; text-decoration: none; }
        .pdf-info a:hover { text-decoration: underline; }
        .occurrences-toggle { cursor: pointer; color: #007bff; font-weight: bold; margin-top: 5px; display: inline-block; }
        .occurrences-list { display: none; margin-top: 8px; padding-left: 15px; border-left: 2px solid #e9ecef; font-size: 12px; max-height: 150px; overflow-y: auto; }
        .occurrences-list li { margin-bottom: 5px; }
        .no-results { text-align: center; padding: 40px; color: #6c757d; font-size: 1.2em; }
    </style>
</head>
<body>
    <div class="container">
        <h1>Git Repository PDF Search</h1> <!-- Changed heading -->

        <div class="stats">
            <h2>Statistics</h2>
            <p>Refs processed: ${stats['processed_refs']} / ${stats['total_refs']}</p>
            <p>Unique commits processed: ${stats['total_commits']}</p>
            <p>Total PDF instances found: ${stats['total_pdfs']}</p> <!-- Changed label -->
            <p>Unique PDFs found: ${stats['unique_pdfs']}</p> <!-- Changed label -->
        </div>

        <div class="controls">
            <div class="search">
                <label for="searchInput">Search by filename or path:</label>
                <input type="text" id="searchInput" placeholder="e.g., report.pdf or documents/archive">
            </div>
            <div class="filter-grid">
                <div class="filter-group">
                    <label for="branchFilter">Branch/Ref:</label>
                    <select id="branchFilter">
                        <option value="all">All Branches/Refs</option>
                    </select>
                </div>
                <div class="filter-group">
                    <label for="sortBy">Sort by:</label>
                    <select id="sortBy">
                        <option value="path">File Path (A-Z)</option>
                        <option value="path-desc">File Path (Z-A)</option>
                        <option value="date-desc">Last Seen (Newest First)</option>
                        <option value="date-asc">First Seen (Oldest First)</option>
                        <option value="occurrences-desc">Occurrences (Most First)</option>
                        <option value="occurrences-asc">Occurrences (Fewest First)</option>
                        <option value="size-desc">Size (Largest First)</option>
                        <option value="size-asc">Size (Smallest First)</option>
                    </select>
                </div>
            </div>
        </div>

        <div class="pdf-grid" id="pdfGrid"> <!-- Changed grid ID -->
''')

        all_branches_refs = set()
        pdf_data_for_js = [] # Renamed variable

        for file_hash, metadata in unique_pdfs.items():
            if not metadata['occurrences']: continue

            occurrences = sorted(metadata['occurrences'], key=lambda x: x['date'], reverse=True)
            latest_occurrence = occurrences[0]
            earliest_occurrence = occurrences[-1]

            latest_path = Path(latest_occurrence['path'])
            latest_branch_dir_name = latest_occurrence['branch'].replace('/', '_').replace('\\', '_')
            latest_commit_str = f"{datetime.datetime.fromisoformat(latest_occurrence['date']).strftime('%Y%m%d_%H%M%S')}_{latest_occurrence['commit_short']}"

            hashed_filename = f"{file_hash}{latest_path.suffix}"
            # Changed base directory in relative path
            relative_pdf_path = Path('pdfs') / latest_branch_dir_name / latest_commit_str / hashed_filename
            display_pdf_src = relative_pdf_path.as_posix()

            # Removed dimensions_str
            size_kb = metadata['size_bytes'] / 1024

            for occ in occurrences:
                all_branches_refs.add(occ['branch'])

            pdf_item_data = { # Renamed variable
                'hash': file_hash,
                'path': metadata['original_path'],
                'branches': list(set(occ['branch'] for occ in occurrences)),
                'occurrences_count': len(occurrences),
                'first_seen': earliest_occurrence['date'],
                'last_seen': latest_occurrence['date'],
                'size_bytes': metadata['size_bytes']
            }
            pdf_data_for_js.append(pdf_item_data) # Renamed variable

            # --- Add Thumbnail HTML ---
            thumbnail_html = ""
            # Check if 'thumbnail_path' exists and is not None in the metadata
            if metadata.get('thumbnail_path'):
                thumbnail_html = f'<img src="{metadata["thumbnail_path"]}" alt="PDF Thumbnail" style="max-width: 100%; height: 150px; object-fit: contain; margin-bottom: 10px;">'
            else:
                 # Placeholder if no thumbnail
                 thumbnail_html = '<div style="height: 150px; display: flex; align-items: center; justify-content: center; background: #f0f0f0; color: #aaa; margin-bottom: 10px; font-size: 12px; border-bottom: 1px solid #eee;">No Thumbnail</div>'

            # --- Write PDF Item HTML ---
            f.write(f'''
        <div class="pdf-item"
             data-hash="{file_hash}"
             data-path="{metadata['original_path'].lower()}"
             data-branches='{json.dumps(list(set(occ['branch'] for occ in occurrences)))}'
             data-occurrences="{len(occurrences)}"
             data-first-seen="{earliest_occurrence['date']}"
             data-last-seen="{latest_occurrence['date']}"
             data-size="{metadata['size_bytes']}">
            {thumbnail_html} <!-- Insert thumbnail or placeholder -->
            <div class="pdf-link-area">
                <a href="{display_pdf_src}" target="_blank" title="Click to open PDF in new tab">
                    Open: {latest_path.name}
                </a>
            </div>
            <div class="pdf-info"> <!-- Changed class name -->
                <p><strong>File:</strong> {latest_path.name}</p>
                <p><strong>Path:</strong> <span class="path">{latest_path.parent}</span></p>
                <p><strong>Size:</strong> {size_kb:.1f} KB</p>
                <!-- Removed Dimensions -->
                <p><strong>Unique Occurrences:</strong> {len(occurrences)}</p>
                <p><strong>First Seen:</strong> {datetime.datetime.fromisoformat(earliest_occurrence['date']).strftime('%Y-%m-%d')}</p>
                <p><strong>Last Seen:</strong> {datetime.datetime.fromisoformat(latest_occurrence['date']).strftime('%Y-%m-%d')}</p>
                <details>
                    <summary class="occurrences-toggle">Show History ({len(occurrences)})</summary>
                    <ul class="occurrences-list">''')

            for occ in occurrences:
                 occ_date_str = datetime.datetime.fromisoformat(occ['date']).strftime('%Y-%m-%d %H:%M')
                 occ_branch_dir = occ['branch'].replace('/', '_').replace('\\', '_')
                 occ_commit_str = f"{datetime.datetime.fromisoformat(occ['date']).strftime('%Y%m%d_%H%M%S')}_{occ['commit_short']}"
                 occ_hashed_filename = f"{file_hash}{Path(occ['path']).suffix}"
                 # Changed base directory in path
                 occ_pdf_path = Path('pdfs') / occ_branch_dir / occ_commit_str / occ_hashed_filename
                 f.write(f'<li><a href="{occ_pdf_path.as_posix()}" target="_blank">{occ_date_str}</a> ({occ["commit_short"]})<br>Ref: {occ["ref"]}<br>Path: {occ["path"]}</li>')

            f.write('''
                    </ul>
                </details>
            </div>
        </div>''')

        f.write('''
        </div>
        <div id="noResults" class="no-results" style="display: none;">No PDFs match your criteria.</div> <!-- Changed text -->
    </div>

    <script>
        const pdfGrid = document.getElementById('pdfGrid'); // Changed grid ID
        const pdfItems = Array.from(pdfGrid.querySelectorAll('.pdf-item')); // Changed item selector
        const searchInput = document.getElementById('searchInput');
        const branchFilter = document.getElementById('branchFilter');
        const sortBySelect = document.getElementById('sortBy');
        const noResultsDiv = document.getElementById('noResults');

        // Populate branch filter
        const branches = ''')
        f.write(json.dumps(sorted(list(all_branches_refs))))
        f.write(''';
        branches.forEach(branch => {
            const option = document.createElement('option');
            option.value = branch;
            option.textContent = branch;
            branchFilter.appendChild(option);
        });

        function updateDisplay() {
            const searchValue = searchInput.value.toLowerCase().trim();
            const branchValue = branchFilter.value;
            const sortValue = sortBySelect.value;

            // Filter PDFs
            let visiblePdfs = pdfItems.filter(item => { // Renamed variable
                const itemPath = item.dataset.path;
                const itemBranches = JSON.parse(item.dataset.branches);

                const matchesSearch = !searchValue || itemPath.includes(searchValue);
                const matchesBranch = branchValue === 'all' || itemBranches.includes(branchValue);

                const isVisible = matchesSearch && matchesBranch;
                item.style.display = isVisible ? 'flex' : 'none'; // Use flex for pdf-item display
                return isVisible;
            });

            // Sort visible PDFs
            visiblePdfs.sort((a, b) => { // Renamed variable
                const pathA = a.dataset.path;
                const pathB = b.dataset.path;
                const lastSeenA = new Date(a.dataset.lastSeen);
                const lastSeenB = new Date(b.dataset.lastSeen);
                const firstSeenA = new Date(a.dataset.firstSeen);
                const firstSeenB = new Date(b.dataset.firstSeen);
                const occurrencesA = parseInt(a.dataset.occurrences);
                const occurrencesB = parseInt(b.dataset.occurrences);
                const sizeA = parseInt(a.dataset.size);
                const sizeB = parseInt(b.dataset.size);

                switch(sortValue) {
                    case 'path': return pathA.localeCompare(pathB);
                    case 'path-desc': return pathB.localeCompare(pathA);
                    case 'date-desc': return lastSeenB - lastSeenA;
                    case 'date-asc': return firstSeenA - firstSeenB;
                    case 'occurrences-desc': return occurrencesB - occurrencesA;
                    case 'occurrences-asc': return occurrencesA - occurrencesB;
                    case 'size-desc': return sizeB - sizeA;
                    case 'size-asc': return sizeA - sizeB;
                    default: return 0;
                }
            });

            // Reorder DOM elements
            visiblePdfs.forEach(item => { // Renamed variable
                pdfGrid.appendChild(item);
            });

            // Show/hide 'no results' message
            noResultsDiv.style.display = visiblePdfs.length === 0 ? 'block' : 'none'; // Renamed variable
        }

        function debounce(func, wait) {
            let timeout;
            return function executedFunction(...args) {
                const later = () => {
                    clearTimeout(timeout);
                    func(...args);
                };
                clearTimeout(timeout);
                timeout = setTimeout(later, wait);
            };
        }

        searchInput.addEventListener('input', debounce(updateDisplay, 300));
        branchFilter.addEventListener('change', updateDisplay);
        sortBySelect.addEventListener('change', updateDisplay);

        updateDisplay();
    </script>
</body>
</html>''')

if __name__ == "__main__":
    import argparse

    default_repo_path = '.'
    # Changed default output dir
    default_output_dir = './git_pdf_results'

    # Changed description
    parser = argparse.ArgumentParser(description='Extract PDFs from a Git repository history.')
    parser.add_argument('repo_path', nargs='?', default=default_repo_path,
                        help=f'Path to the local Git repository (default: {default_repo_path})')
    # Changed output help text and default
    parser.add_argument('--output', default=default_output_dir,
                        help=f'Directory to store results (default: {default_output_dir})')

    args = parser.parse_args()

    repo_abs_path = Path(args.repo_path).resolve()
    output_abs_path = Path(args.output).resolve()

    print(f"Repository to analyze: {repo_abs_path}")
    print(f"Output directory: {output_abs_path}")

    # Changed function call
    extract_git_pdfs(str(repo_abs_path), str(output_abs_path))

    ## to run:
    # python ImageExtraction/git_pdf_search.py --output ImageExtraction/git_pdf_results
    # start ImageExtraction/git_pdf_results/index.html
