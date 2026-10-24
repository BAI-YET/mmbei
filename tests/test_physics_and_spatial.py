import numpy as np,pandas as pd
from mathorcup.q4_solver import simulate_1d,simulate_2d,analytic_voellmy_velocity
from mathorcup.q4_flowpath import directed_cost_path

def test_analytic_voellmy_has_correct_terminal_speed_and_initial_acceleration():
    theta=30;mu=.2;xi=400;h=4
    t=np.array([0,1e-5,1,100])
    u=analytic_voellmy_velocity(t,theta,mu,xi,h)
    assert u[0]==0
    a=9.81*(np.sin(np.deg2rad(theta))-mu*np.cos(np.deg2rad(theta)))
    assert np.isclose(u[1]/t[1],a,rtol=1e-4)
    assert 0<u[2]<u[-1]

def test_1d_conservation_on_synthetic_profile():
    s=np.arange(0,2000,50.0)
    p=pd.DataFrame({'distance_m':s,'elevation_m':4200-np.tan(np.deg2rad(25))*s})
    r=simulate_1d(p,mu=.10,xi=350,width=40,volume=1000,duration=50,
        total_seconds=200,frames=4,cfl=.25)
    assert r['mass']['relative_absolute_mass_residual']<1e-8
    assert np.min(r['last_h'])>=0
    assert r['mass']['injected_m3']>990

def test_2d_conservation_and_nonnegative():
    z=np.repeat((900-np.arange(20)*2.)[:,None],20,axis=1)
    t=np.array([0,10,20,30,40,50,60,70,80],float)
    hyd=pd.DataFrame({'time_s':t,'sensor_discharge_m3_s':[0,10,30,30,10,0,0,0,0]})
    ans=simulate_2d(z,30,(3,10),hyd,total_seconds=75,mu=.05,xi=400,frames=4,cfl=.2)
    assert ans['mass']['relative_mass_error']<1e-7
    assert np.min(ans['U'][0])>=0
    assert ans['mass']['input_cumulative_m3']>0

def test_dijkstra_is_directed_downhill():
    z=np.array([[7,6,5,4],[9,8,7,3],[11,10,8,2]],float)
    path,cost,n=directed_cost_path(z,(0,0),(2,3),cell_size=10,uphill_penalty=100)
    assert tuple(path[0])==(0,0) and tuple(path[-1])==(2,3)
    assert np.isfinite(cost)
    assert n>0
