config.title = "simple_injector"
print(config.title)
config.dimensions = 2
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
M_inj = 2.0 --start with sonic injection (choking occurs at injector)
--freestream flow conditions
Pinf = 7.69e3 --landsberg
Patm=Pinf/20
inf_T=254
inf_dummy = FlowState:new{p=Pinf, T=inf_T,massf=massf_air} --landsberg
inf_vx = M_inf * inf_dummy.a
inf = FlowState:new{p=Pinf, T=inf_T, velx = inf_vx, massf=massf_air} -- freestream flow condition at combustor

initial = FlowState:new{p=Patm, T=inf_T, velx=0,vely=0,massf=massf_air} --initial state with zero velocity to be set upon each block
--Injector conditions
--inj_P = 1.0e+7
inj_Pstag=1611.0e3
inj_Tstag= 300.0
gamma = 1.4

inj_T = inj_Tstag / (1+0.2*M_inj*M_inj)
inj_P = inj_Pstag*math.pow(inj_T/inj_Tstag,3.5)

inj_dummy = FlowState:new{p=inj_P, T=inj_T, massf=massf_fuel} --used storage pressure and tmp of hydrogen gas
inj_vy = M_inj * inj_dummy.a
inj = FlowState:new{p=inj_P, T=inj_T, vely = inj_vy, massf=massf_fuel}
--inj = FlowState:new{p=inj_P, T=inj_T, vely = 0}
print("Vx at the inlet is	: ",inf_vx)
print("Vy at the injector is	: ",inj_vy)


-- Demo: Verify Mach number of inflow and compute dynamic pressure.

--print("inflow=", inflow)
--print("T=", inf.T, "density=", inf.rho, "sound speed= ", inf.a, "Mach = ", inf.M)
--print("inflow Mach number=", 1000.0/inflow.a)
--print("dynamic pressure q=", 1/2*inflow.rho*1.0e6)
-- Set up two quadrilaterals in the (x,y)-plane by first defining
-- the corner nodes, then the lines between those corners.
--x and y units are in m
ymax=20.0e-3
xmax=100.0e-3
xinj=30.0e-3
injw = 0.2e-3
a0 = Vector3:new{x=0.0, y=0.0}
a1 = Vector3:new{x=xinj, y=0.0}
a2 = Vector3:new{x=(xinj+injw), y=0.}
a3 = Vector3:new{x=xmax, y=0.}
b0 = Vector3:new{x=0.0, y=ymax}
b1 = Vector3:new{x=xinj, y=ymax}
b2 = Vector3:new{x=(xinj+injw), y=ymax}
b3 = Vector3:new{x=xmax, y=ymax}

--lower boundary
a0a1 = Line:new{p0=a0, p1=a1} 
a1a2 = Line:new{p0=a1, p1=a2} 
a2a3 = Line:new{p0=a2, p1=a3}
--upper boundary
b0b1 = Line:new{p0=b0, p1=b1} 
b1b2 = Line:new{p0=b1, p1=b2} 
b2b3 = Line:new{p0=b2, p1=b3}
--vertical lines
a0b0 = Line:new{p0=a0, p1=b0} 
a1b1 = Line:new{p0=a1, p1=b1} 
a2b2 = Line:new{p0=a2, p1=b2}
a3b3 = Line:new{p0=a3, p1=b3}

----new small block underneath which its north boundary acts as injector

quad0 = CoonsPatch:new{north=b0b1, east=a1b1, south=a0a1, west=a0b0}
quad1 = CoonsPatch:new{north=b1b2, east=a2b2, south=a1a2, west=a1b1}
quad2 = CoonsPatch:new{north=b2b3, east=a3b3, south=a2a3, west=a2b2}
--do some clustering
rcfR = RobertsFunction:new{end0=false, end1=true, beta=1.02}
rcfL = RobertsFunction:new{end0=true, end1=false, beta=1.02}
rcfD = RobertsFunction:new{end0=true, end1=false, beta=1.01}
cfList1={east=rcfD,west=rcfD} -- clusters blk1
cfList0={north=rcfR,east=rcfD,south=rcfR,west=rcfD}-- clusters blk0
cfList2={north=rcfL,east=rcfD,south=rcfL,west=rcfD} -- clusters blk2

-- Mesh the patches, with particular discretisation.
nx0 = 30; nx1 = 10; nx2 = 600; ny = 125
blk0 = StructuredGrid:new{psurface=quad0, niv=nx0+1, njv=ny+1,cfList=cfList0}
blk1 = StructuredGrid:new{psurface=quad1, niv=nx1+1, njv=ny+1,cfList=cfList1}
blk2 = StructuredGrid:new{psurface=quad2, niv=nx2+1, njv=ny+1,cfList=cfList2}
-- Define the flow-solution blocks.
F_inf = FluidBlock:new{grid=blk0, initialState= initial}
F_inj = FluidBlock:new{grid=blk1, initialState= initial}
F_out = FluidBlock:new{grid=blk2, initialState = initial}
identifyBlockConnections()
-- Set boundary conditions.
F_inf.bcList[west] = InFlowBC_Supersonic:new{flowState=inf}
--F_inf.bcList[east] = OutFlowBC_Simple:new{}
--F_inj.bcList[west] = InFlowBC_ConstFlux:new{flowState=inf}
--F_inf.bcList[south] = InFlowBC_ShockFitting:new{flowState=inj}
F_inf.bcList[south] = WAllBC_NoSlip_Adiabatic:new{}
F_inj.bcList[south] = InFlowBC_Supersonic:new{flowState=inj}
--F_inj.bcList[south] = InFlowBC_FromStagnation:new{stagnationState=inj}
--F_inj.bcList[north] = WallBC_WithSlip:new{}
--F_out.bcList[west] = InFlowBC_Supersonic:new{flowState=inf}
F_out.bcList[east] = OutFlowBC_Simple:new{}
F_out.bcList[south]=WallBC_NoSlip_Adiabatic:new{}
--history points on the cell with error ([0.0301, 0.0114586, 0])
setHistoryPoint{x=0.0301,y=0.0114586,z=0}

max_t=0.03
max_step=10000000
-- Do a little more setting of global data.
config.flux_calculator="adaptive_efm_ausmdv"
config.gasdynamic_update_scheme='euler'
config.cfl_value = 0.5
config.moving_grid=false
--config.viscous=true
--config.turbulence_model="adaptive"
config.thermo_interpolator="pT"
config.max_invalid_cells=20
config.adjust_invalid_cell_data=true
config.max_time = max_t -- seconds
config.max_step = max_step
config.dt_init = 1.0e-9
config.dt_plot = 1e-5
--config.dt_history = config.max_time/100.0
config.dt_history = 1.0e-8
--config.dt_loads = 1.0e-8 --calculate loads at evert 1.0e-5 seconds
--dofile("sketch-domain.lua")
