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
M_inj = 2.0  --start with sonic injection (choking occurs at injector)
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
l	=1.5e-3
ymax	=30.0e-3
xmax	=100.0e-3
xinj	=30.00e-3
zmax=ymax
print("ready for grid")
x0=0
x1=xinj
x2=xinj+l
x3=xinj+l+0.5*(xmax-xinj-l)
x4=xmax

y1=0.5*l
y2=(ymax-0.5*l)/2
y3=ymax

aa00 = Vector3:new{x=x0, 	y=0.0,	z=0.0}
aa01 = Vector3:new{x=x1, 	y=0.0,	z=0.0}
aa02 = Vector3:new{x=x2, 	y=0.0,	z=0.0}
aa03 = Vector3:new{x=x3, 	y=0.0,	z=0.0}
aa04 = Vector3:new{x=x4, 	y=0.0,	z=0.0}
ab00 = Vector3:new{x=x0, 	y=y1,	z=0.0}
ab01 = Vector3:new{x=x1, 	y=y1,	z=0.0}
ab02 = Vector3:new{x=x2, 	y=y1,	z=0.0}
ab03 = Vector3:new{x=x3, 	y=y1,	z=0.0}
ab04 = Vector3:new{x=x4, 	y=y1,	z=0.0}
ac00 = Vector3:new{x=x0, 	y=y2,	z=0.0}
ac01 = Vector3:new{x=x1, 	y=y2,	z=0.0}
ac02 = Vector3:new{x=x2, 	y=y2,	z=0.0}
ac03 = Vector3:new{x=x3, 	y=y2,	z=0.0}
ac04 = Vector3:new{x=x4, 	y=y2,	z=0.0}
ad00 = Vector3:new{x=x0, 	y=y3,	z=0.0}
ad01 = Vector3:new{x=x1, 	y=y3,	z=0.0}
ad02 = Vector3:new{x=x2, 	y=y3,	z=0.0}
ad03 = Vector3:new{x=x3, 	y=y3,	z=0.0}
ad04 = Vector3:new{x=x4, 	y=y3,	z=0.0}

ba00 = Vector3:new{x=x0, 	y=0.0,	z=zmax/2}
ba01 = Vector3:new{x=x1, 	y=0.0,	z=zmax/2}
ba02 = Vector3:new{x=x2, 	y=0.0,	z=zmax/2}
ba03 = Vector3:new{x=x3, 	y=0.0,	z=zmax/2}
ba04 = Vector3:new{x=x4, 	y=0.0,	z=zmax/2}
bb00 = Vector3:new{x=x0, 	y=y1,	z=zmax/2}
bb01 = Vector3:new{x=x1, 	y=y1,	z=zmax/2}
bb02 = Vector3:new{x=x2, 	y=y1,	z=zmax/2}
bb03 = Vector3:new{x=x3, 	y=y1,	z=zmax/2}
bb04 = Vector3:new{x=x4, 	y=y1,	z=zmax/2}
bc00 = Vector3:new{x=x0, 	y=y2,	z=zmax/2}
bc01 = Vector3:new{x=x1, 	y=y2,	z=zmax/2}
bc02 = Vector3:new{x=x2, 	y=y2,	z=zmax/2}
bc03 = Vector3:new{x=x3, 	y=y2,	z=zmax/2}
bc04 = Vector3:new{x=x4, 	y=y2,	z=zmax/2}
bd00 = Vector3:new{x=x0, 	y=y3,	z=zmax/2}
bd01 = Vector3:new{x=x1, 	y=y3,	z=zmax/2}
bd02 = Vector3:new{x=x2, 	y=y3,	z=zmax/2}
bd03 = Vector3:new{x=x3, 	y=y3,	z=zmax/2}
bd04 = Vector3:new{x=x4, 	y=y3,	z=zmax/2}

ca00 = Vector3:new{x=x0, 	y=0.0,	z=zmax}
ca01 = Vector3:new{x=x1, 	y=0.0,	z=zmax}
ca02 = Vector3:new{x=x2, 	y=0.0,	z=zmax}
ca03 = Vector3:new{x=x3, 	y=0.0,	z=zmax}
ca04 = Vector3:new{x=x4, 	y=0.0,	z=zmax}
cb00 = Vector3:new{x=x0, 	y=y1,	z=zmax}
cb01 = Vector3:new{x=x1, 	y=y1,	z=zmax}
cb02 = Vector3:new{x=x2, 	y=y1,	z=zmax}
cb03 = Vector3:new{x=x3, 	y=y1,	z=zmax}
cb04 = Vector3:new{x=x4, 	y=y1,	z=zmax}
cc00 = Vector3:new{x=x0, 	y=y2,	z=zmax}
cc01 = Vector3:new{x=x1, 	y=y2,	z=zmax}
cc02 = Vector3:new{x=x2, 	y=y2,	z=zmax}
cc03 = Vector3:new{x=x3, 	y=y2,	z=zmax}
cc04 = Vector3:new{x=x4, 	y=y2,	z=zmax}
cd00 = Vector3:new{x=x0, 	y=y3,	z=zmax}
cd01 = Vector3:new{x=x1, 	y=y3,	z=zmax}
cd02 = Vector3:new{x=x2, 	y=y3,	z=zmax}
cd03 = Vector3:new{x=x3, 	y=y3,	z=zmax}
cd04 = Vector3:new{x=x4, 	y=y3,	z=zmax}
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
nx0 = 30; nx1 = 20; nx2=60; ny0 = 10; ny1 = 30;  nz = 30

v0	= TFIVolume:new{vertices={aa00,aa01,ab01,ab00,ba00,ba01,bb01,bb00}}
v1	= TFIVolume:new{vertices={aa01,aa02,ab02,ab01,ba01,ba02,bb02,bb01}}
v2	= TFIVolume:new{vertices={aa02,aa03,ab03,ab02,ba02,ba03,bb03,bb02}}
v3	= TFIVolume:new{vertices={aa03,aa04,ab04,ab03,ba03,ba04,bb04,bb03}}
v4	= TFIVolume:new{vertices={ab00,ab01,ac01,ac00,bb00,bb01,bc01,bc00}}
v5	= TFIVolume:new{vertices={ab01,ab02,ac02,ac01,bb01,bb02,bc02,bc01}}
v6	= TFIVolume:new{vertices={ab02,ab03,ac03,ac02,bb02,bb03,bc03,bc02}}
v7	= TFIVolume:new{vertices={ab03,ab04,ac04,ac03,bb03,bb04,bc04,bc03}}
v8	= TFIVolume:new{vertices={ac00,ac01,ad01,ad00,bc00,bc01,bd01,bd00}}
v9	= TFIVolume:new{vertices={ac01,ac02,ad02,ad01,bc01,bc02,bd02,bd01}}
v10	= TFIVolume:new{vertices={ac02,ac03,ad03,ad02,bc02,bc03,bd03,bd02}}
v11	= TFIVolume:new{vertices={ac03,ac04,ad04,ad03,bc03,bc04,bd04,bd03}}

v12	= TFIVolume:new{vertices={ba00,ba01,bb01,bb00,ca00,ca01,cb01,cb00}}
v13	= TFIVolume:new{vertices={ba01,ba02,bb02,bb01,ca01,ca02,cb02,cb01}}
v14	= TFIVolume:new{vertices={ba02,ba03,bb03,bb02,ca02,ca03,cb03,cb02}}
v15	= TFIVolume:new{vertices={ba03,ba04,bb04,bb03,ca03,ca04,cb04,cb03}}
v16	= TFIVolume:new{vertices={bb00,bb01,bc01,bc00,cb00,cb01,cc01,cc00}}
v17	= TFIVolume:new{vertices={bb01,bb02,bc02,bc01,cb01,cb02,cc02,cc01}}
v18	= TFIVolume:new{vertices={bb02,bb03,bc03,bc02,cb02,cb03,cc03,cc02}}
v19	= TFIVolume:new{vertices={bb03,bb04,bc04,bc03,cb03,cb04,cc04,cc03}}
v20	= TFIVolume:new{vertices={bc00,bc01,bd01,bd00,cc00,cc01,cd01,cd00}}
v21	= TFIVolume:new{vertices={bc01,bc02,bd02,bd01,cc01,cc02,cd02,cd01}}
v22	= TFIVolume:new{vertices={bc02,bc03,bd03,bd02,cc02,cc03,cd03,cd02}}
v23	= TFIVolume:new{vertices={bc03,bc04,bd04,bd03,cc03,cc04,cd04,cd03}}
print("produced volumes")
blk0 = StructuredGrid:new{pvolume=v0, niv=nx0+1, njv=ny0+1,nkv=nz+1,
cfList={edge01=cf1_1,edge12=nil,edge32=cf1_1,edge03=nil,edge45=cf1_1,edge56=nil,edge76=cf1_1,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk1 = StructuredGrid:new{pvolume=v1, niv=nx1+1, njv=ny0+1,nkv=nz+1,
cfList={edge01=nil,edge12=nil,edge32=nil,edge03=nil,edge45=nil,edge56=nil,edge76=nil,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk2 = StructuredGrid:new{pvolume=v2, niv=nx2+1, njv=ny0+1,nkv=nz+1,
cfList={edge01=cf0_1,edge12=nil,edge32=cf0_1,edge03=nil,edge45=cf0_1,edge56=nil,edge76=cf0_1,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk3 = StructuredGrid:new{pvolume=v3, niv=nx0+1, njv=ny1+1,nkv=nz+1,
cfList={edge01=nil,edge12=nil,edge32=nil,edge03=nil,edge45=nil,edge56=nil,edge76=nil,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
print("produced 0-3")
blk4 = StructuredGrid:new{pvolume=v4, niv=nx1+1, njv=ny1+1,nkv=nz+1,
cfList={edge01=cf1_1,edge12=nil,edge32=cf1_3,edge03=nil,edge45=cf1_1,edge56=nil,edge76=cf1_3,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk5 = StructuredGrid:new{pvolume=v5, niv=nx2+1, njv=ny1+1,nkv=nz+1,
cfList={edge01=nil,edge12=nil,edge32=nil,edge03=nil,edge45=nil,edge56=nil,edge76=nil,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk6 = StructuredGrid:new{pvolume=v6, niv=nx0+1, njv=ny0+1,nkv=nz+1,
cfList={edge01=cf1_1,edge12=nil,edge32=cf1_1,edge03=nil,edge45=cf1_1,edge56=nil,edge76=cf1_1,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk7 = StructuredGrid:new{pvolume=v7, niv=nx1+1, njv=ny0+1,nkv=nz+1,
cfList={edge01=nil,edge12=nil,edge32=nil,edge03=nil,edge45=nil,edge56=nil,edge76=nil,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
print("produced 0-7")
blk8 = StructuredGrid:new{pvolume=v8, niv=nx2+1, njv=ny0+1,nkv=nz+1,
cfList={edge01=cf1_3,edge12=nil,edge32=cf1_5,edge03=nil,edge45=cf1_3,edge56=nil,edge76=cf1_5,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk9 = StructuredGrid:new{pvolume=v9, niv=nx0+1, njv=ny1+1,nkv=nz+1,
cfList={edge01=nil,edge12=nil,edge32=nil,edge03=nil,edge45=nil,edge56=nil,edge76=nil,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk10 = StructuredGrid:new{pvolume=v10, niv=nx1+1, njv=ny1+1,nkv=nz+1,
cfList={edge01=cf0_1,edge12=nil,edge32=cf0_1,edge03=nil,edge45=cf0_1,edge56=cf0_3,edge76=nil,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk11 = StructuredGrid:new{pvolume=v11, niv=nx2+1, njv=ny1+1,nkv=nz+1,
cfList={edge01=nil,edge12=nil,edge32=nil,edge03=nil,edge45=nil,edge56=nil,edge76=nil,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
print("produced 0-11")
blk12 = StructuredGrid:new{pvolume=v12, niv=nx0+1, njv=ny0+1,nkv=nz+1,
cfList={edge01=cf1_1,edge12=nil,edge32=cf1_1,edge03=nil,edge45=cf1_1,edge56=nil,edge76=cf1_1,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk13 = StructuredGrid:new{pvolume=v13, niv=nx1+1, njv=ny0+1,nkv=nz+1,
cfList={edge01=nil,edge12=nil,edge32=nil,edge03=nil,edge45=nil,edge56=nil,edge76=nil,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk14 = StructuredGrid:new{pvolume=v2, niv=nx2+1, njv=ny0+1,nkv=nz+1,
cfList={edge01=cf0_1,edge12=nil,edge32=cf0_1,edge03=nil,edge45=cf0_3,edge56=nil,edge76=cf0_3,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk15 = StructuredGrid:new{pvolume=v15, niv=nx0+1, njv=ny1+1,nkv=nz+1,
cfList={edge01=nil,edge12=nil,edge32=nil,edge03=nil,edge45=nil,edge56=nil,edge76=nil,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
print("produced 0-15")
blk16 = StructuredGrid:new{pvolume=v16, niv=nx1+1, njv=ny1+1,nkv=nz+1,
cfList={edge01=cf1_1,edge12=nil,edge32=cf1_3,edge03=nil,edge45=cf1_1,edge56=nil,edge76=cf1_3,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk17 = StructuredGrid:new{pvolume=v17, niv=nx2+1, njv=ny1+1,nkv=nz+1,
cfList={edge01=nil,edge12=nil,edge32=nil,edge03=nil,edge45=nil,edge56=nil,edge76=nil,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk18 = StructuredGrid:new{pvolume=v18, niv=nx0+1, njv=ny0+1,nkv=nz+1,
cfList={edge01=cf0_1,edge12=nil,edge32=cf0_1,edge03=nil,edge45=cf0_3,edge56=nil,edge76=cf0_3,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk19 = StructuredGrid:new{pvolume=v19, niv=nx1+1, njv=ny0+1,nkv=nz+1,
cfList={edge01=nil,edge12=nil,edge32=nil,edge03=nil,edge45=nil,edge56=nil,edge76=nil,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
print("produced 0-19")
blk20 = StructuredGrid:new{pvolume=v20, niv=nx2+1, njv=ny0+1,nkv=nz+1,
cfList={edge01=cf1_3,edge12=nil,edge32=cf1_5,edge03=nil,edge45=cf1_3,edge56=nil,edge76=cf1_5,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk21 = StructuredGrid:new{pvolume=v21, niv=nx0+1, njv=ny1+1,nkv=nz+1,
cfList={edge01=nil,edge12=nil,edge32=nil,edge03=nil,edge45=nil,edge56=nil,edge76=nil,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk22 = StructuredGrid:new{pvolume=v22, niv=nx1+1, njv=ny1+1,nkv=nz+1,
cfList={edge01=cf0_1,edge12=nil,edge32=cf0_3,edge03=nil,edge45=cf0_3,edge56=nil,edge76=cf0_3,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
blk23 = StructuredGrid:new{pvolume=v23, niv=nx2+1, njv=ny1+1,nkv=nz+1,
cfList={edge01=nil,edge12=nil,edge32=nil,edge03=nil,edge45=nil,edge56=nil,edge76=nil,edge47=nil,edge04=rcfz,edge15=rcfz,edge26=rcfz,edge37=rcfz}}
identifyBlockConnections()
print("produced 0-23")
-- Define the flow-solution blocks. Naming convention: F_(flow type)_blocknumber(0~23)
--inlet blocks
F_inf_0 = FluidBlock:new{grid=blk0, initialState= initial,
bcList={west=InFlowBC_Supersonic:new{flowState=inf}}}
F_inf_4 = FluidBlock:new{grid=blk4, initialState= initial,
bcList={west=InFlowBC_Supersonic:new{flowState=inf}}}
F_inf_8 = FluidBlock:new{grid=blk8, initialState= initial,
bcList={west=InFlowBC_Supersonic:new{flowState=inf}}}
F_inf_12 = FluidBlock:new{grid=blk12, initialState= initial,
bcList={west=InFlowBC_Supersonic:new{flowState=inf}}}
F_inf_16 = FluidBlock:new{grid=blk16, initialState= initial,
bcList={west=InFlowBC_Supersonic:new{flowState=inf}}}
F_inf_20 = FluidBlock:new{grid=blk20, initialState= initial,
bcList={west=InFlowBC_Supersonic:new{flowState=inf}}}

F_inj_1 = FluidBlock:new{grid=blk1, initialState=initial,
bcList={bottom=InFlowBC_Supersonic:new{flowState=inj}}}

F_inf_2 = FluidBlock:new{grid=blk2, initialState= initial}
F_inf_5 = FluidBlock:new{grid=blk5, initialState= initial}
F_inf_6 = FluidBlock:new{grid=blk6, initialState= initial}
F_inf_9 = FluidBlock:new{grid=blk9, initialState= initial}
F_inf_10 = FluidBlock:new{grid=blk10, initialState= initial}
F_inf_13 = FluidBlock:new{grid=blk13, initialState= initial}
F_inf_14 = FluidBlock:new{grid=blk14, initialState= initial}
F_inf_17 = FluidBlock:new{grid=blk17, initialState= initial}
F_inf_18 = FluidBlock:new{grid=blk18, initialState= initial}
F_inf_21 = FluidBlock:new{grid=blk21, initialState= initial}
F_inf_22 = FluidBlock:new{grid=blk22, initialState= initial}

F_out_3 = FluidBlock:new{grid=blk3, initialState= initial,
bcList={west=OutFlowBC_Simple:new{}}}
F_out_7 = FluidBlock:new{grid=blk7, initialState= initial,
bcList={west=OutFlowBC_Simple:new{}}}
F_out_11 = FluidBlock:new{grid=blk11, initialState= initial,
bcList={west=OutFlowBC_Simple:new{}}}
F_out_15 = FluidBlock:new{grid=blk15, initialState= initial,
bcList={west=OutFlowBC_Simple:new{}}}
F_out_19 = FluidBlock:new{grid=blk19, initialState= initial,
bcList={west=OutFlowBC_Simple:new{}}}
F_out_23 = FluidBlock:new{grid=blk23, initialState= initial,
bcList={west=OutFlowBC_Simple:new{}}}

--F_side= FluidBlock:new{grid=blk4, initialState= initial}
--F_out = FluidBlock:new{grid=blk2, initialState = initial}
--F_out2 = FluidBlock:new{grid=blk5, initialState = initial}
print("so far so good...")
identifyBlockConnections()
-- Set boundary conditions.

--F_inf.bcList[south] = InFlowBC_Supersonic:new{flowState=inf}
--F_inf2.bcList[south] = InFlowBC_Supersonic:new{flowState=inf}
--F_inj.bcList[bottom] = InFlowBC_Supersonic:new{flowState=inj}
--F_out.bcList[north] = OutFlowBC_Simple:new{}
--F_out2.bcList[north] = OutFlowBC_Simple:new{}
max_t=0.001
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
--config.max_step = max_step
config.dt_init = 1.0e-7
config.dt_plot = 1e-5
--config.dt_history = config.max_time/100.0
--config.dt_history = 1.0e-8
--config.dt_loads = 1.0e-8 --calculate loads at evert 1.0e-5 seconds
