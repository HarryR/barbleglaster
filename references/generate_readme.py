import bibtexparser
import os

def generate_readme(bibtex_file="ref.bib", output_file="README.md"):
    bib_database = bibtexparser.parse_file(bibtex_file)
    readme_content = "# References\n\nThis directory contains references cited in this project.\n\n"
    for entry in bib_database.entries:
        readme_content += f"## {entry['title']}\n\n"
        readme_content += f"- **type**: {entry.entry_type}\n"

        # Add all properties
        for key, value in sorted(entry.items()):
            # Skip ENTRYTYPE as we've already added it
            if key == 'ENTRYTYPE' or key == 'ID':
                continue

            if key == 'file':
                # Extract filename from file field and create link
                file_path = value.strip('{}').split(':')[0]
                file_name = os.path.basename(file_path)
                readme_content += f"- **file**: [{file_name}]({file_path})\n"
            elif key == 'url':
                readme_content += f"- **url**: [{value}]({value})\n"
            elif key == 'biburl':
                readme_content += f"- **biburl**: [{value}]({value})\n"
            elif key == 'doi':
                # Check if DOI is already a URL, if not, make it one
                if value.startswith('http'):
                    doi_url = value
                else:
                    doi_url = f"https://doi.org/{value}"
                readme_content += f"- **doi**: [{value}]({doi_url})\n"
            else:
                readme_content += f"- **{key}**: {value}\n"

        readme_content += "\n"

    # Write to README.md
    with open(output_file, 'w', encoding='utf-8') as readme_file:
        readme_file.write(readme_content)

if __name__ == "__main__":
    generate_readme()
