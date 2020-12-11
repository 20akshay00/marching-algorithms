inf = open("points.vtk")
stripped_lines = [l.lstrip() for l in inf.readlines()]
inf.close()

# write the new, stripped lines to a file
outf = open("points.vtk", "w")
outf.write("".join(stripped_lines))
outf.close()
