import sys
import subprocess
from pathlib import Path


def run():
    if len(sys.argv) < 2:
        print("Usage: python run_pipeline.py <A>")
        sys.exit(1)

    A = sys.argv[1]
    
    ROOT = Path(__file__).resolve().parent.parent
    BIN = ROOT / "scripts"/"build"
    SCRIPTS = ROOT / "scripts"
    INPUTS = ROOT / "inputs"
    LIB = ROOT / "cgshop26" / "XOR-SAT"

    # 1) convert_format.py A
    print("Step 1: convert_format.py")
    subprocess.run(
        ["python3", SCRIPTS / "convert_format.py", A],
        check=True
    )

    # 2) 2sat A
    print("Step 2: 2sat")
    subprocess.run(
        [BIN / "2sat_linux", A],
        check=True
    )

    # 3) merge.py
    print("Step 3: merge.py")
    subprocess.run(
        ["python3", SCRIPTS / "merge.py"],
        check=True
    )

    # 4) cryptominisat5
    print("Step 4: cryptominisat5")
    with open(SCRIPTS / "output.txt", "w") as out:
        result = subprocess.run(
            [
                LIB / "cryptominisat5",
                "--verb", "0",
                "--threads", "16",
                ROOT / "total.cnf"
            ],
            stdout=out
        )

    sol = "s"
    if result.returncode == 10:
        print("SAT")
        s = "SAT"
    elif result.returncode == 20:
        print("UNSAT")
        s = "UNSAT"
    else:
        s = "ERROR"
        print("Solver error, exit code:", result.returncode)
        sys.exit(1)

    with open("verdict.txt", "w") as f:
        f.write(str(s))

    print("Pipeline finished successfully")


if __name__ == "__main__":
    run()

