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

    m = int(sys.argv[2])

    init_dist = 5
    if len(sys.argv) >= 4:
        init_dist = int(sys.argv[3])

    matrix = [[0 for _ in range(m)] for _ in range(m)]
    distances = [0 for _ in range(m)]
    lista_i = [2,1,3,5,0,4]

    init_dist = [[0 for _ in range(m)] for _ in range(m)]
    init_dist[0][1] = 38
    init_dist[0][2] = 34
    init_dist[0][3] = 39
    init_dist[0][4] = 45
    init_dist[0][5] = 45
    init_dist[1][2] = 18
    init_dist[1][3] = 23
    init_dist[1][4] = 46
    init_dist[1][5] = 35
    init_dist[2][3] = 19
    init_dist[2][4] = 42
    init_dist[2][5] = 31
    init_dist[3][4] = 45
    init_dist[3][5] = 36
    init_dist[4][5] = 46

    for i in range(m):
        for j in range(i):
             init_dist[i][j] = init_dist[j][i]


    for i in lista_i:
        for j in range(m):
            if matrix[i][j] != 0:
                continue
            sol = [-1,-1]
            cnt = 0
            sacuvaj_listu_u_fajl([2,i,j], "../config/" + "/selection.txt")
            l = 0
            r = init_dist[i][j]

            while l <= r:
                cnt += 1
                if cnt > 6:
                    break
                mid = (l+r)//2
                mid = int(mid)
                d1 = mid//2
                d1 = int(d1)
                d2 = d1
                d2 = int(d2)
                if mid%2==1:
                    d2 += 1
                    d2 = int(d2)

                distances[i] = d1
                distances[j] = d2
                sacuvaj_listu_u_fajl(distances, "../config/distances.txt")

                print(str(i) + " " + str(j) + " " + str(d1) + " " + str(d2)+ " " + str(mid))
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
                        sol = [d1, d2]
                        r = mid-1
                    else:
                        l = mid+1
                

            d1 = sol[0]
            d2 = sol[1]
            matrix[i][j] = d1 + d2
            matrix[j][i] = d1 + d2
            with open("./" + str(i) + "_" + str(j) + ".txt", "w") as f:
                f.write(str(matrix[i][j]))

    with open("./" + A + "_pairwise.txt", "w") as f:
        for row in matrix:
            f.write(" ".join(map(str, row)) + "\n")




if __name__ == "__main__":
    run()
