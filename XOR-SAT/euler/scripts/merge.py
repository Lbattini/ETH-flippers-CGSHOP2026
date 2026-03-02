from pathlib import Path

def append_cnf_files(folder_path, output_file):
    folder = Path(folder_path)

    # Get all .cnf files and sort numerically by filename
    cnf_files = sorted(
        folder.glob("*.cnf"),
        key=lambda p: int(p.stem)
    )

    with open(output_file, "w", encoding="utf-8") as out:
        for cnf_file in cnf_files:
            with open(cnf_file, "r", encoding="utf-8") as f:
                out.write(f.read())

                # Optional: add a newline between files
                out.write("\n")

if __name__ == "__main__":
    append_cnf_files("../formulations", "../total.cnf")
