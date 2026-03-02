import json
import os
import sys


def convert(json_path):
    # izlazni fajl sa istim imenom, ali .in ekstenzija
    base = os.path.splitext(json_path)[0]
    out_path = base + ".in"

    with open(json_path, "r") as f:
        data = json.load(f)

    points_x = data["points_x"]
    points_y = data["points_y"]
    triangulations = data["triangulations"]

    lines = []

    # placeholder za prvi red (popunjavamo kasnije)
    lines.append("")

    # --- Points: x_i y_i ---
    for x, y in zip(points_x, points_y):
        lines.append(f"{x} {y}")

    lines.append("")  # prazan red

    ed_cnt = 0
    # --- Triangulations: grane kao indeksi (u v) ---
    for tri in triangulations:
        cnt = 0
        for edge in tri:
            cnt += 1
            u, v = edge
            lines.append(f"{u} {v}")
        ed_cnt = cnt  # pretpostavka: sve triangulacije imaju isti broj grana
        lines.append("")

    tri_cnt = len(triangulations)

    # prvi red: broj tacaka, broj grana, broj triangulacija
    lines[0] = f"{len(points_x)} {ed_cnt} {tri_cnt}"

    # --- upis u .in fajl ---
    with open(out_path, "w") as f:
        f.write("\n".join(lines))

    print(f"Generated: {out_path}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python convert.py <json_file>")
    else:
        convert("../data/" + sys.argv[1] + ".json")
