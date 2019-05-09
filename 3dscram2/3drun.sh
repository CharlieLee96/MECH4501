e4shared --prep --job=3d_grid3_2
e4shared --run --job=3d_grid3_2
e4shared --post --job=3d_grid3_2 --vtk-xml --add-vars="mach" --tindx-plot=all
e4shared --custom-post --script-file=locate-bow-shock.lua

