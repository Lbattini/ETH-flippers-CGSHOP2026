import sys
import os
import subprocess
from pathlib import Path

def ucitaj_fajl_u_listu(putanja):
    rezultat = []

    with open(putanja, "r") as f:
        for red in f:
            brojevi = red.strip().split()
            
            brojevi = [int(x) for x in brojevi]
            rezultat.append(brojevi)

    return rezultat

def sacuvaj_listu_u_fajl(lista, putanja):
    with open(putanja, "w") as f:
        f.write(" ".join(str(x) for x in lista))

def procitaj_i_smanji(putanja):
    with open(putanja, "r") as f:
        broj = int(f.read().strip())

    if broj < 0:
        sys.exit(0)

    # upis smanjenog broja
    with open(putanja, "w") as f:
        f.write(str(broj - 1))

    return broj



def run():
    if len(sys.argv) < 2:
        print("Usage: python run_pipeline.py <argument A>")
        sys.exit(1)

    A = sys.argv[1]

    l = ucitaj_fajl_u_listu("../distances.txt")


    for lista in l:
    	sacuvaj_listu_u_fajl(lista, "../config/" + "/distances.txt")
    	sacuvaj_listu_u_fajl([0], "../config/" + "/selection.txt")

    	print("Step 0: run_all.py")
    	subprocess.run(
        	[
            	"python3",
            	r"../scripts/run_all.py",
            	A
        	],
        	check=True
    	)

    	with open("verdict.txt", "r") as f:
            s = f.read().strip().split()
            if s[0] == "SAT":
                print("FINISHED, SOLUTION EXISTS")
                exit(0)
if __name__ == "__main__":
    run()
