config.title = "sphere in argon flow, M = 6.0"
config.dimensions = 2
config.axisymmetric = true

-- Flow conditions
--setGasModel("ideal-argon-gas.lua")
--setGasModel("ideal-air.lua")
-- Given Mach number, but need velocity (m/s)

setGasModel('H2air-model.lua') 	-- fuel

massf_air={N2=0.78,O2=0.22} 	--mass flow rate of air
massf_fuel={H2=1.0}		--mass flow rate of hydrogen gas

p_inf = 10.0e3 -- Pa
T_inf = 300.0 -- K
-- Use temporary flow state to have Eilmer compute sound speed, a
temp = FlowState:new{p=p_inf, T=T_inf}
a_inf = temp.a
M_inf = 6.0
u_inf = M_inf*a_inf


p_int=101000
T_int=298
-- Set flow conditions
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf,massf=massf_air}
--initial = FlowSolution:new{jobName='sphere', dir="../coarse",
                           --tindx='last', nBlocks=1}
initial=FlowState:new{p=p_int, T=T_int, velx=0}

M_f=1
P_f=500.0e3
T_f= 295.0
inj_dummy = FlowState:new{p=P_f, T=T_f} --used storage pressure and tmp of hydrogen gas
Vy_f = M_f * inj_dummy.a

fuel = FlowState:new{p=p_f, T=T_f, vely=Vy_f , massf=massf_fuel}
print(Vy_f)
-- Specify geometry
r = 10.0e-3 -- m
A = Vector3:new{x=0.0, y=0.0}
B = Vector3:new{x=r,   y=r  }
C = Vector3:new{x=r,   y=0.0}
D = Vector3:new{x=-0.3*r, y=0.0}
E = Vector3:new{x=r, y=2.0*r}
d0 = Vector3:new{x=D.x, y=0.7*r}
d1 = Vector3:new{x=D.x + 0.7*(E.x - D.x), y=E.y}

-- Set boundary paths
DA = Line:new{p0=D, p1=A}
AB = Arc:new{p0=A, p1=B, centre=C}
EB = Line:new{p0=E, p1=B}
DE = Bezier:new{points={D, d0, d1, E}}

-- Build patch, grid and block
factor = 2
niCells = factor*30
njCells = factor*40
quad0 = makePatch{north=EB, east=AB, south=DA, west=DE}
grid0 = StructuredGrid:new{psurface=quad0, niv=niCells+1, njv=njCells+1}
blk0 = FluidBlockArray{grid=grid0, initialState=initial, nib=1, njb=4,
                      bcList={west=InFlowBC_Supersonic:new{flowState=inflow},
                              north=OutFlowBC_Simple:new{},
                              south=InFlowBC_Supersonic:new{flowState=fuel}}
			      --east =WallBC_NoSlip:new{}}
}

-- Set some simulation parameters
config.flux_calculator = "ausmdv"
--config.max_time = 3.0*r/u_inf
config.max_time = 8e-5
config.max_step = 500000
config.dt_init = 1.0e-8
config.cfl_value = 0.5
config.dt_plot = config.max_time/10.0

-- Set some locations for recording detailed data
setHistoryPoint{x=A.x, y=A.y}
setHistoryPoint{x=B.x, y=B.y}
config.dt_history = 1.0e-7 
config.viscous=true

-- Set eilmer to record surface loads
config.compute_loads = true
config.dt_loads = 1.0e-5














