config.title = "simple_injector"
print(config.title)
config.dimensions = 2
config.axisymmetric = true
setGasModel('H2air-model.lua') 	-- fuel

massf_air={N2=0.78,O2=0.22} 	--mass flow rate of air
massf_fuel={H2=1.0}		--mass flow rate of hydrogen gas

--inlet conditions are by landsberg.
M_inlet = 4.586 		--combustor freestream flow condition specified by Landsberg at Mach 10
M_injec = 1.0 			--start with sonic injection (choking occurs at injector)
--freestream flow conditions
P_inflow = 50.51e3 --landsberg
P_atm=P_inflow/20
T_free=1100
inf_dummy = FlowState:new{p=P_inflow, T=T_free} --landsberg
Vx_free = M_inlet * inf_dummy.a
inf = FlowState:new{p=P_inflow, T=T_free, velx = Vx_free, massf=massf_air} -- freestream flow condition at combustor

initial = FlowState:new{p=P_atm, T=T_free, velx = 0,vely = 0,massf=massf_air} --initial state with zero velocity to be set upon each block
--Injector conditions
--P_injec = 1.0e+7
Pstag_injec=500.0e3
Tstag_injec= 295.0
gamma = 1.4

T_injec = Tstag_injec / (1+0.2*M_injec*M_injec)
P_injec = Pstag_injec*math.pow(T_injec/Tstag_injec,3.5)

inj_dummy = FlowState:new{p=P_injec, T=T_injec} --used storage pressure and tmp of hydrogen gas
Vy_injec = M_injec * inj_dummy.a
injec = FlowState:new{p=P_injec, T=T_injec, vely = Vy_injec,massf=massf_fuel}
--inj = FlowState:new{p=P_injec, T=T_injec, vely = 0}



-- Demo: Verify Mach number of inflow and compute dynamic pressure.

--print("inflow=", inflow)
--print("T=", inf.T, "density=", inf.rho, "sound speed= ", inf.a, "Mach = ", inf.M)
--print("inflow Mach number=", 1000.0/inflow.a)
--print("dynamic pressure q=", 1/2*inflow.rho*1.0e6)
-- Set up two quadrilaterals in the (x,y)-plane by first defining
-- the corner nodes, then the lines between those corners.
--x and y units are in m
ymax=20.0e-3
xmax=100.0e-3		--maximum x (length of combustor)
xinj=30.0e-3		--x coordinate of start of injector 
injw = 0.2e-3		--width of injector in 2D

a0 = Vector3:new{x=0.0, y=-(injw*4)}			--y is the length of the injector
a1 = Vector3:new{x=xinj, y=-(injw*4)}
a2 = Vector3:new{x=(xinj+injw), y=-(injw*4)}

b0 = Vector3:new{x=0.0, y=0.}
b1 = Vector3:new{x=xinj, y=0.}
b2 = Vector3:new{x=(xinj+injw), y=0.}
b3 = Vector3:new{x=(b2.x)+(1/3)*(xmax-(xinj+injw)), y = 0.}
b4 = Vector3:new{x=xmax, y=0.}

--c0 = Vector3:new{x=0.0, y=ymax}
--c1 = Vector3:new{x=xinj, y=ymax}
--c2 = Vector3:new{x=(xinj+injw), y=ymax}
--c3 = Vector3:new{x=(c2.x)+(1/3)*(xmax-(xinj+injw)), y = ymax}
--c4 = Vector3:new{x=xmax, y=ymax}
c0 = Vector3:new{x=b0.x, y=ymax}
c1 = Vector3:new{x=b1.x, y=ymax}
c2 = Vector3:new{x=b2.x, y=ymax}
c3 = Vector3:new{x=b3.x, y = ymax}
c4 = Vector3:new{x=b4.x, y=ymax}
--injector

a1b1 = Line:new{p0=a1, p1=b1} 
a1a2 = Line:new{p0=a1, p1=a2}
a2b2 = Line:new{p0=a2, p1=b2}

--lower boundary
b0b1 = Line:new{p0=b0, p1=b1} 
b1b2 = Line:new{p0=b1, p1=b2} 
b2b3 = Line:new{p0=b2, p1=b3}
b3b4 = Line:new{p0=b3, p1=b4}
--upper boundary
c0c1 = Line:new{p0=c0, p1=c1} 
c1c2 = Line:new{p0=c1, p1=c2} 
c2c3 = Line:new{p0=c2, p1=c3}
c3c4 = Line:new{p0=c3, p1=c4}
--vertical lines
b0c0 = Line:new{p0=b0, p1=c0} 
b1c1 = Line:new{p0=b1, p1=c1} 
b2c2 = Line:new{p0=b2, p1=c2}
b3c3 = Line:new{p0=b3, p1=c3}
b4c4 = Line:new{p0=b4, p1=c4}
--patches
quad0 = CoonsPatch:new{north=c0c1, east=b1c1, south=b0b1, west=b0c0}
--quad1 = CoonsPatch:new{north=b1b2, east=a2b2, south=a1a2, west=a1b1}
quad2 = CoonsPatch:new{north=c1c2, east=b2c2, south=b1b2, west=b1c1}
quad3 = CoonsPatch:new{north=c2c3, east=b3c3, south=b2b3, west=b2c2}
quad4 = CoonsPatch:new{north=c3c4, east=b4c4, south=b3b4, west=b3c3}

-- Mesh the patches, with particular discretisation.
nx0 = 8; nx2 = 2; nx3 = 20;nx4 = 20 ny = 20; ny_inj = 3
blk0 = StructuredGrid:new{psurface=quad0, niv=nx0+1, njv=ny+1,cfList=cfList0} --cfList0 governs clustering of blk0.
--blk1 = StructuredGrid:new{psurface=quad1, niv=nx2+1, njv=ny_inj+1}
blk2 = StructuredGrid:new{psurface=quad2, niv=nx2+1, njv=ny+1}
blk3 = StructuredGrid:new{psurface=quad2, niv=nx3+1, njv=ny+1,cfList=cfList3} --cfList0 governs clustering of blk3.
blk4 = StructuredGrid:new{psurface=quad2, niv=nx4+1, njv=ny+1,cfList=cfList4}
-- Define the flow-solution blocks. Start with all initial state except F_1 (injector block)
F_0 = FluidBlock:new{grid=blk0, initialState= initial}
--F_1 = FluidBlock:new{grid=blk1, initialState= injec}
F_2 = FluidBlock:new{grid=blk2, initialState = initial}
F_3 = FluidBlock:new{grid=blk3, initialState = initial}
F_4 = FluidBlock:new{grid=blk4, initialState = initial}
identifyBlockConnections()
-- Set boundary conditions.
--blk0
F_0.bcList[west] = InFlowBC_Supersonic:new{flowState=inf}
--blk1 (injector)
--F_1.bcList[north] = OutFlowBC_Simple:new{} -- check wallbc of block 2
--F_2.bcList[south] = InFlowBC_Supersonic:new{flowState=injec}
F_2.bcList[east] = OutFlowBC_Simple:new{}
F_3.bcList[east] = OutFlowBC_Simple:new{}
F_4.bcList[east] = OutFlowBC_Simple:new{}
--do some clustering
rcfR = RobertsFunction:new{end0=false, end1=true, beta=1.08}
rcfL = RobertsFunction:new{end0=true, end1=false, beta=1.05}
--cfList0={east=rcfL,south=rcfR}		-- clusters blk0
cfList3={west=rcfL,south=rcfL} 		-- clusters blk2
cfList3={west=rcfL,south=rcfL}
max_t=0.001
max_step=100000
-- Do a little more setting of global data.
max_invalid_cells = 10
config.flux_calculator="adaptive_efm_ausmdv"
config.gasdynamic_update_scheme='euler'
config.cfl_value = 5
config.moving_grid=false
config.viscous=true
config.turbulence_model="adaptive"
config.max_time = max_t -- seconds
config.max_step = max_step
config.dt_init = 1.0e-13
--config.dt_init = 1.0e-6
config.dt_plot = 1e-3
--config.dt_history = config.max_time/100.0
--config.dt_history = 10.0e-7
--dofile("sketch-domain.lua")
