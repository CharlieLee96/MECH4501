config.title = "3d_injector"
print(config.title)
config.dimensions = 3
config.axisymmetric = false
-- The gas model is defined via a gas-model file.
--nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
--print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)

setGasModel('H2air-model.lua')
massf_air={N2=0.78,O2=0.22}
massf_fuel={H2=1.0}
----setting up flow conditions
--Mach numbers
M_inf = 5.63 --combustor freestream flow condition specified by Landsberg at Mach 10
M_inj = 2 --start with sonic injection (choking occurs at injector)
--freestream flow conditions
Pinf = 7.69e3 --landsberg
Patm=Pinf/20
inf_T=254
inf_dummy = FlowState:new{p=Pinf, T=inf_T,massf=massf_air} --landsberg
inf_vx = M_inf * inf_dummy.a
inf = FlowState:new{p=Pinf, T=inf_T, vely = inf_vx, massf=massf_air} -- freestream flow condition at combustor
print("T at inlet	=",inf.T)
print("P at inlet	=",inf.p)
print("rho at inlet	=",inf.rho)
initial = FlowState:new{p=Patm, T=inf_T, velx=0,velz=0,massf=massf_air} --initial state with zero velocity to be set upon each block
print("initial_done") 
--Injector conditions
--inj_P = 1.0e+7
inj_Pstag=1611.0e3
inj_Tstag= 300.0
gamma = 1.4

inj_T = inj_Tstag / (1+0.2*M_inj*M_inj)
inj_P = inj_Pstag*math.pow(inj_T/inj_Tstag,3.5)
print("Tstag at injector=",inj_T)
print("Pstag at injector=",inj_P)
inj_dummy = FlowState:new{p=inj_P, T=inj_T, massf=massf_fuel} --used storage pressure and tmp of hydrogen gas
inj_vz = M_inj * inj_dummy.a
inj = FlowState:new{p=inj_P, T=inj_T, velz = inj_vz, massf=massf_fuel}
print("T at injector	=",inj.T)
print("P at injector	=",inj.p)
print("rho at injector	=",inj.rho)
--inj = FlowState:new{p=inj_P, T=inj_T, velz = 0}
print("Vx at the inlet is	: ",inf_vx)
print("Vz at the injector is	: ",inj_vz)
--w	=5.542e-3
w	=1.639e-3
ymax	=20.0e-3
xmax	=100.0e-3
xinj	=30.00e-3
yinj	=w
injw 	= 4*w
H=ymax
print("ready for grid")

aa00 = Vector3:new{x=0.0, 	y=0.0,	z=0.0}
aa01 = Vector3:new{x=xinj, 	y=0.0,	z=0.0}
aa02 = Vector3:new{x=(xinj+4*w), y=0.0,	z=0.0}
aa03 = Vector3:new{x=xmax, 	y=0.0,	z=0.0}
ab00 = Vector3:new{x=0.0, 	y=yinj,z=0.0}
ab01 = Vector3:new{x=xinj, 	y=yinj,z=0.0}
ab02 = Vector3:new{x=(xinj+4*w), y=yinj,z=0.0}
ab03 = Vector3:new{x=xmax, 	y=yinj,z=0.0}
ac00 = Vector3:new{x=0.0, 	y=ymax,z=0.0}
ac01 = Vector3:new{x=xinj, 	y=ymax,z=0.0}
ac02 = Vector3:new{x=(xinj+4*w), y=ymax,z=0.0}
ac03 = Vector3:new{x=xmax, 	y=ymax,z=0.0}

ba00 = Vector3:new{x=0.0, 	y=0.0,	z=H}
ba01 = Vector3:new{x=xinj, 	y=0.0,	z=H}
ba02 = Vector3:new{x=(xinj+4*w), y=0.0,	z=H}
ba03 = Vector3:new{x=xmax, 	y=0.0,	z=H}
bb00 = Vector3:new{x=0.0, 	y=yinj,z=H}
bb01 = Vector3:new{x=xinj, 	y=yinj,z=H}
bb02 = Vector3:new{x=(xinj+4*w), y=yinj,z=H}
bb03 = Vector3:new{x=xmax, 	y=yinj,z=H}
bc00 = Vector3:new{x=0.0, 	y=ymax,z=H}
bc01 = Vector3:new{x=xinj, 	y=ymax,z=H}
bc02 = Vector3:new{x=(xinj+4*w), y=ymax,z=H}
bc03 = Vector3:new{x=xmax, 	y=ymax,z=H}
--define clustering a bit (1 is strongest)

cf0_1 = RobertsFunction:new{end0=true,end1=false,beta=1.1}
cf0_3 = RobertsFunction:new{end0=true,end1=false,beta=1.3}
cf0_5 = RobertsFunction:new{end0=true,end1=false,beta=1.5}
cf1_1 = RobertsFunction:new{end0=false,end1=true,beta=1.1}
cf1_3 = RobertsFunction:new{end0=false,end1=true,beta=1.3}
cf1_5 = RobertsFunction:new{end0=false,end1=true,beta=1.5}

rcfL = RobertsFunction:new{end0=true, end1=false, beta=1.02}
rcfR = RobertsFunction:new{end0=false, end1=true, beta=1.02}
rcfz = RobertsFunction:new{end0=true, end1=false, beta=1.01}


--3d blocks
nx0 = 30; nx1 = 20; nx2=60; ny0 = 10; ny1 = 30;  nz = 60

v0l={aa00,aa01,ab01,ab00,ba00,ba01,bb01,bb00}
v1l={aa01,aa02,ab02,ab01,ba01,ba02,bb02,bb01}
v2l={aa02,aa03,ab03,ab02,ba02,ba03,bb03,bb02}
v3l={ab00,ab01,ac01,ac00,bb00,bb01,bc01,bc00}
v4l={ab01,ab02,ac02,ac01,bb01,bb02,bc02,bc01}
v5l={ab02,ab03,ac03,ac02,bb02,bb03,bc03,bc02}
print(v0l)
print(v1l)
print(v2l)
print(v3l)
print(v4l)
print(v5l)
v0	= TFIVolume:new{vertices=v0l}
v1	= TFIVolume:new{vertices=v1l}
v2	= TFIVolume:new{vertices=v2l}
v3	= TFIVolume:new{vertices=v3l}
v4	= TFIVolume:new{vertices=v4l}
v5	= TFIVolume:new{vertices=v5l}


blk0 = StructuredGrid:new{pvolume=v0, niv=nx0+1, njv=ny0+1,nkv=nz+1,
cfList={edge01=cf1_1,edge12=nil,edge32=cf1_1,edge03=nil,edge45=cf1_1,edge56=nil,edge76=cf1_1,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk1 = StructuredGrid:new{pvolume=v1, niv=nx1+1, njv=ny0+1,nkv=nz+1,
cfList={edge01=nil,edge12=nil,edge32=nil,edge03=nil,edge45=nil,edge56=nil,edge76=nil,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk2 = StructuredGrid:new{pvolume=v2, niv=nx2+1, njv=ny0+1,nkv=nz+1,
cfList={edge01=cf0_1,edge12=nil,edge32=cf0_1,edge03=nil,edge45=cf0_1,edge56=nil,edge76=cf0_1,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}

blk3 = StructuredGrid:new{pvolume=v3, niv=nx0+1, njv=ny1+1,nkv=nz+1,
cfList={edge01=cf1_1,edge12=cf0_3,edge32=cf1_1,edge03=nil,edge45=cf1_1,edge56=nil,edge76=cf1_1,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk4 = StructuredGrid:new{pvolume=v4, niv=nx1+1, njv=ny1+1,nkv=nz+1,
cfList={edge01=nil,edge12=cf0_3,edge32=nil,edge03=cf0_3,edge45=nil,edge56=cf0_5,edge76=nil,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk5 = StructuredGrid:new{pvolume=v5, niv=nx2+1, njv=ny1+1,nkv=nz+1,
cfList={edge01=cf0_1,edge12=nil,edge32=cf0_1,edge03=cf0_3,edge45=cf0_1,edge56=nil,edge76=cf0_1,edge47=cf0_5,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}


-- Define the flow-solution blocks.
F_inf = FluidBlockArray{grid=blk0, initialState= initial,bcList={west=InFlowBC_Supersonic:new{flowState=inf}},nib=1,njb=1,nkb=1}
F_inf2 = FluidBlockArray{grid=blk3, initialState= initial,bcList={west=InFlowBC_Supersonic:new{flowState=inf}},nib=1,njb=3,nkb=1}
F_inj = FluidBlockArray{grid=blk1, initialState= initial,bcList={bottom=InFlowBC_Supersonic:new{flowState=inj}},nib=1,njb=1,nkb=1}
F_side= FluidBlockArray{grid=blk4, initialState= initial,nib=1,njb=3,nkb=1}
F_out = FluidBlockArray{grid=blk2, initialState = initial,bcList={east=OutFlowBC_Simple:new{}},nib=2,njb=1,nkb=1}
F_out2 = FluidBlockArray{grid=blk5, initialState = initial,bcList={east=OutFlowBC_Simple:new{}},nib=2,njb=3,nkb=1}
print("so far so good...")
identifyBlockConnections()
-- Set boundary conditions.
--F_inf.bcList[south] = InFlowBC_Supersonic:new{flowState=inf}
--F_inf2.bcList[south] = InFlowBC_Supersonic:new{flowState=inf}
--F_inj.bcList[bottom] = InFlowBC_Supersonic:new{flowState=inj}
--F_out.bcList[north] = OutFlowBC_Simple:new{}
--F_out2.bcList[north] = OutFlowBC_Simple:new{}
max_t=0.001
max_step=10000000
-- Do a little more setting of global data.
config.flux_calculator="adaptive_efm_ausmdv"
config.gasdynamic_update_scheme='euler'
config.cfl_value = 0.25
config.moving_grid=false
config.viscous=true
--config.turbulence_model="adaptive"
config.thermo_interpolator="pT"
config.max_invalid_cells=200
config.adjust_invalid_cell_data=false
config.max_time = max_t -- seconds
config.max_step = max_step
config.dt_init = 1.0e-7
config.dt_plot = 1e-5
--config.dt_history = config.max_time/100.0
--config.dt_history = 1.0e-8
--config.dt_loads = 1.0e-8 --calculate loads at evert 1.0e-5 seconds
