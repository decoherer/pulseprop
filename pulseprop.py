import numpy as np
from numpy import pi,exp,log,abs,sqrt,conj
from wavedata import Wave,wrange,logrange,timeit,track
import scipy, scipy.integrate
from sellmeier import index,polingperiod,qpmwavelengths,groupvelocity,groupindex,groupvelocitydispersion,dispersionlength
# from joblib import Memory
# memory = Memory('c:/backup', verbose=0) # use as @memory.cache

def transformlimitedbandwidth(dt,λ): # dt=tFWHM in ns, returns Δλ in nm
    df = 2*log(2)/pi/dt # in GHz
    return df*λ**2/299792458
def transformlimitedpulseamplitude(f,f0,dt): # tFWHM * fFWHM = 2 log(2) / pi = 0.441271
    return exp(-pi**2 * (f-f0)**2 * dt**2 / log(2))
def timebandwidthproduct(A,λ0,fwhm=False): # returns 1 for transform-limited gaussian pulse
    B = temporal2spectral(A,λ0,vsfreq=1)
    return pi/2/log(2) * A.magsqr().fwhm() * B.magsqr().fwhm() if fwhm else 4*pi*sqrt( A.magsqr().secondmoment() * B.magsqr().secondmoment() )
def peakpower(J,τ): # J = pulse energy in nJ, τ = tFWHM of gaussian in ns
    return J / τ /sqrt(pi/log(16)) # peak power in W
def pulsesum(x,t,f0,reptime,dt,sell='air'): # dt = tFWHM
    c = 299.792458 # in mm/ns
    num,frr = 10*int(reptime/dt),reptime
    fs = f0 + frr*np.linspace(-num,num,2*num+1)
    ns = (lambda x:1)(fs)
    return transformlimitedpulseamplitude(fs,f0,dt) * exp(1j*2*pi*fs*(x*ns/c-t))
def coupledwaveshg(η,L,dk,P1,P2=0,tol=1e-6): # η in %/W/cm², L in mm, input power P1 in mW, dk in µm⁻¹
    g = sqrt(η/100/1000/10**2) # g in units of 1/√mW/mm
    E0s = [sqrt(P1)+0j,sqrt(P2)+0j] # initial fields, units of √mW (Eᵢ≡√Pᵢ)
    def ODEs(z,Es):
        E1,E2 = Es
        return [-1j * g * exp(-1j*dk*1000*z) * E2 * conj(E1),
                -1j * g * exp(+1j*dk*1000*z) * E1 * E1 ]
    sol = scipy.integrate.solve_ivp(ODEs, [0,L], E0s, dense_output=True, rtol=tol, atol=tol)
    z = np.linspace(0,L,2001)
    return [Wave(e,z) for e in sol.sol(z)]
def coupledwavesfg(η,L,dk,λ1,λ2,E1=0,E2=0,E3=0,loss1=0,loss2=0,loss3=0,tol=1e-6):
    # η in %/W/cm², L in mm, input fields in √mW, dk in µm⁻¹, loss in dB/cm
    λ3 = 1/(1/λ1+1/λ2)
    h = sqrt(η/100/1000/10**2) # h in units of 1/√mW/mm
    E0s = [E+0j for E in (E1,E2,E3)] # initial fields, units of √mW (Eᵢ≡√Pᵢ)
    α1,α2,α3 = [0.01*loss*np.log(10) for loss in (loss1,loss2,loss3)] # e**-α = 10**(-0.1*loss) → α = 0.1*loss*ln10, extra 0.1 for dB/mm
    def ODEs(z,Es):
        E1,E2,E3 = Es
        return [-1j * h * exp(-1j*dk*1000*z) * E3 * conj(E2) * λ3/λ1 - 0.5*α1*E1,
                -1j * h * exp(-1j*dk*1000*z) * E3 * conj(E1) * λ3/λ2 - 0.5*α2*E2,
                -1j * h * exp(+1j*dk*1000*z) * E1 * E2  - 0.5*α3*E3]
    sol = scipy.integrate.solve_ivp(ODEs, [0,L], E0s, dense_output=True, rtol=tol, atol=tol)
    z = np.linspace(0,L,2001)
    return [Wave(e,z) for e in sol.sol(z)]
def singlepulseprop(τ,L,P,λ0=1550,sell='ktp',dt=0.005,plot=False): # τ = fwhm pulse duration in ps, L in mm, P = peak power in W
    c = 0.299792458 # mm/ps = m/ns
    v1,v3 = c*groupvelocity(λ0,sell),c*groupvelocity(λ0/2,sell)
    print(v1,v3)
    ts = wrange(-20*τ,int(2*L/c)+20*τ,dt) # ps 
    P1 = Wave(P*exp(-4*log(2)*ts**2/τ**2),ts)
    P3 = P1 + 0
    E1 = sqrt(P1) + 0j
    E3 = sqrt(P3) + 0j
    def propagate(A,x,λ):
        ωs = 2*pi*np.fft.fftfreq(len(A), dt) # ωs in rad/ps
        ω0 = 1e6*2*pi*c/λ
        λs = 1e6*2*pi*c/(ω0+ωs) # in nm
        op = -1j*index(λs,sell)*(ω0+ωs)/c + 1j*index(λ,sell)*ω0/c
        def ode(z,A):
            return np.fft.ifft(op*np.fft.fft(A))
        sol = scipy.integrate.solve_ivp(ode, [0,x], A)
        return Wave(sol.y[:,-1], A.x)
    A1 = propagate(E1,L,λ0)
    A3 = propagate(E3,L,λ0/2)
    if plot:
        Wave.plots(P1,P3,A1.magsqr().rename(1),A3.magsqr().rename(3),lines=[0,L/v1,L/v3],x='time (ps)',y='power (W)',xlim=(L/v1-2*τ,L/v3+2*τ),ylim=(-0.1,1.1),grid=1,
            save=f'pulse propagation {τ:g}ps {L:g}mm {P:g}W {λ0:g}nm {sell}')
    return A1,A3
def pulseenergy(A): # A = amplitude in √W vs time in ps, λ0 in nm # returns pulse energy in pJ
    return A.magsqr().sum()*A.dx()
def temporal2spectral(A,λ0,vsfreq=False,correctphaseramp=False,plot=False,logplot=True): # A = amplitude in √W vs time in ps, λ0 in nm
    A = A if not correctphaseramp else A.removefftphaseramp() # correct for linear phase due to carrier frequency
    if plot:
        xlim = 'f' if logplot else (A.magsqr().maxloc()-2*A.magsqr().fwhm(),A.magsqr().maxloc()+2*A.magsqr().fwhm())
        A.magsqr().plot(x='time (ps)',y='power (W)',seed=λ0,xlim=xlim,legendtext=f'{A.magsqr().sum()*A.dx():g}pJ pulse energy',grid=1,log=logplot)
    dt,N = A.dx(),len(A)
    df,fs = 1/dt/N, np.fft.fftshift(np.fft.fftfreq(N,dt)) # in THz # print('Δf (GHz)',1000*(fs[1]-fs[0]),'\nΔf (GHz)',1000*df)
    c = 0.299792458  # mm/ps = m/ns
    f0 = 1e6*c/λ0 # in THz
    λs,dλ0 = 1e6*c/(f0+fs), 1e-6*c*df/f0**2 # in nm
    B = Wave(dt*np.fft.fftshift(np.fft.fft(A)),fs) # y in √(J·ps), x in THz
    B = B if vsfreq else Wave(1e3*sqrt(c)*B.y/λs,λs) # y in √J/nm, x in nm
    if plot:
        xlim = 'f' if logplot else (B.magsqr().maxloc()-2*B.magsqr().fwhm(),B.magsqr().maxloc()+2*B.magsqr().fwhm())
        if vsfreq: # print(f"{B.magsqr().area():g}pJ")
            B.magsqr().plot(x='frequency (THz)',y='energy (pJ/THz)',seed=λ0,xlim=xlim,legendtext=f'{B.magsqr().area():g}pJ pulse energy',grid=1,log=logplot)
        else:
            B.magsqr().plot(x='wavelength (nm)',y='power (pJ/nm)',seed=λ0,xlim=xlim,legendtext=f'{B.magsqr().sum()*dλ0*1e12:g}pJ pulse energy',grid=1,log=logplot)
    return B
# @timeit
# @memory.cache
def shgpulseprop(τ,L,P,η=0,λ0=1550,sell='ktp',Type='zzz',Λ=None,dt=0.005,d=None,rtol=1e-3,atol=1e-6,plot=False): # τ = fwhm pulse duration in ps, L in mm, P = peak power in W, η in %/W/cm², Λ in µm
    c = 0.299792458  # mm/ps = m/ns
    λ1,λ3 = λ0,λ0/2
    v1,v3 = c*groupvelocity(λ1,sell),c*groupvelocity(λ3,sell) # print(v1,v3)
    kg = 0 if Λ is None else 1e3*2*pi * (1/polingperiod(λ1,λ1,sell,Type) - 1/Λ) # in mm⁻¹
    d = d if d is not None else lambda z:1
    g = sqrt(0.01*0.01*η) # 1/100x for %, 1/100x for cm² to mm²
    ts = wrange(-20*τ,int(L/v3/dt)*dt+20*τ,dt)  # ps
    P1 = Wave(P*exp(-4*log(2)*ts**2/τ**2),ts)
    P3 = 0*P1 + 0
    N,E1,E3 = len(P1),sqrt(P1) + 0j,sqrt(P3) + 0j
    ωs = 2*pi*np.fft.fftfreq(N, dt)  # ωs in rad/ps
    ω1,ω3 = 1e6*2*pi*c/λ1,1e6*2*pi*c/λ3
    λ1s,λ3s = 1e6*2*pi*c/(ω1+ωs),1e6*2*pi*c/(ω3+ωs)  # in nm
    prop1 = -1j*index(λ1s,sell)*(ω1+ωs)/c + 1j*index(λ1,sell)*ω1/c
    prop3 = -1j*index(λ3s,sell)*(ω3+ωs)/c + 1j*index(λ3,sell)*ω3/c
    def ode(z, y):
        E1z,E3z = y[:N],y[N:]
        dE1 = np.fft.ifft(prop1*np.fft.fft(E1z))
        dE3 = np.fft.ifft(prop3*np.fft.fft(E3z))
        dE1 += 1j * g * d(z/L) * E3z * np.conj(E1z) * exp(1j*kg*z)
        dE3 += 1j * g * conj(d(z/L)) * E1z**2 * exp(-1j*kg*z)
        return np.concatenate([dE1, dE3])
    sol = scipy.integrate.solve_ivp(ode, [0,L], np.concatenate([E1.y, E3.y]), atol=atol, rtol=rtol, t_eval=[0,L], method='DOP853')
    A1,A3 = Wave(sol.y[:N,-1], ts),Wave(sol.y[N:,-1], ts)
    # print(f" Pin {P1.sum()*P1.dx():g}pJ, Pout {A1.magsqr().sum()*A1.dx():g}pJ, Pshg {A3.magsqr().sum()*A3.dx():g}pJ")
    if plot:
        Wave.plots(P1,A1.magsqr(),
                lines=[0,L/v1,L/v3],x='time (ps)',y='pump peak power (W)',grid=1,log=0,
                xlim=(L/v1-2*τ,L/v3+2*τ),ylim=(-0.1*P,1.1*P),
                rightwaves=[A3.magsqr()],rightlim=(-0.1*A3.magsqr().max(),1.1*A3.magsqr().max()),rightlabel='SHG peak power (W)',
                save=f'pulse propagation {τ:g}ps {L:g}mm {P:g}W {λ0:g}nm {sell}')
        Wave.plots(P1,A1.magsqr().rename('pump'),A3.magsqr().rename('SHG'),
                lines=[0,L/v1,L/v3],x='time (ps)',y='power (W)',grid=1,log=1,
                xlim=(L/v1-2*τ,L/v3+2*τ),ylim=(1e-4,1e5),
                save=f'pulse propagation {τ:g}ps {L:g}mm {P:g}W {λ0:g}nm {sell} log')
        Bin,Bout,Bshg = temporal2spectral(E1,λ0),temporal2spectral(A1,λ0),temporal2spectral(A3,λ0/2)
        xlim = (Bin.maxloc()-2*Bin.fwhm(),Bin.maxloc()+2*Bin.fwhm())
        Wave.plots(Bin.magsqr().rename('pump in'),Bout.magsqr().rename('pump out'),x='λ (nm)',xlim=xlim,save=f'shg pump pulse spectrum {τ:g}ps {L:g}mm {P:g}W {λ0:g}nm {sell}')
        Wave.plots(Wave(Bin.y,Bin.x/2).magsqr().rename('pump in'),Bshg.magsqr().rename('SHG'),x='λ (nm)',xlim=(xlim[0]/2,xlim[1]/2),save=f'shg output pulse spectrum {τ:g}ps {L:g}mm {P:g}W {λ0:g}nm {sell}')
    return A1,A3
# @timeit
# @memory.cache
def sfgpulseprop(τ, L, P1, P2, η=0, λ1=1550, λ2=775, sell='ktp', Type='zzz', Λ=None, dt=0.005, rtol=1e-3, atol=1e-6, plot=False):
    # τ = fwhm pulse duration in ps, L in mm, P = peak power in W, η in %/W/cm²
    c = 0.299792458 # mm/ps = m/ns
    ω1, ω2 = 1e6*2*np.pi*c/λ1, 1e6*2*np.pi*c/λ2
    λ3,ω3 = 1/(1/λ1 + 1/λ2), ω1 + ω2
    v1, v2, v3 = (c*groupvelocity(λ,sell) for λ in (λ1, λ2, λ3))
    kg = 0 if Λ is None else 1e3*2*pi * (1/polingperiod(λ1,λ2,sell,Type) - 1/Λ) # in mm⁻¹
    g = sqrt(0.01*0.01*η) # 1/100x for %, 1/100x for cm² to mm²
    ts = wrange(-20*τ, int(L/v3/dt)*dt+20*τ, dt)  # ps
    N = len(ts)
    E1 = sqrt(P1 * np.exp(-4*np.log(2)*ts**2/τ**2)) + 0j
    E2 = sqrt(P2 * np.exp(-4*np.log(2)*ts**2/τ**2)) + 0j
    E3 = 0*E1
    ωs = 2*np.pi*np.fft.fftfreq(N, dt)
    λ1s,λ2s,λ3s = [1e6*2*np.pi*c / (ω + ωs) for ω in (ω1,ω2,ω3)]
    prop1 = -1j*index(λ1s,sell)*(ω1+ωs)/c + 1j*index(λ1,sell)*ω1/c
    prop2 = -1j*index(λ2s,sell)*(ω2+ωs)/c + 1j*index(λ2,sell)*ω2/c
    prop3 = -1j*index(λ3s,sell)*(ω3+ωs)/c + 1j*index(λ3,sell)*ω3/c
    def ode(z, y):
        E1z, E2z, E3z = np.split(y, 3)
        dE1 = np.fft.ifft(prop1*np.fft.fft(E1z))
        dE2 = np.fft.ifft(prop2*np.fft.fft(E2z))
        dE3 = np.fft.ifft(prop3*np.fft.fft(E3z))
        dE1 += 1j*g*(λ3/λ1) * E2z.conj()*E3z * np.exp(1j*kg*z)
        dE2 += 1j*g*(λ3/λ2) * E1z.conj()*E3z * np.exp(1j*kg*z)
        dE3 += 1j*g * E1z*E2z * np.exp(-1j*kg*z)
        return np.concatenate([dE1, dE2, dE3])
    y0 = np.concatenate([E1, E2, E3])
    sol = scipy.integrate.solve_ivp(ode, [0,L], y0, atol=atol, rtol=rtol, t_eval=[0,L], method='DOP853') # RK45, BDF, Radau, LSODA, RK23
    # sol = scipy.integrate.solve_ivp(ode, [0,L], y0, atol=atol, rtol=rtol, method='LSODA', t_eval=[0,L], first_step=0.0001, min_step=0.0001, max_step=0.001) # RK45, BDF, Radau, LSODA, RK23
    A1,A2,A3 = (Wave(sol.y[i*N:(i+1)*N, -1], ts) for i in range(3))
    if plot:
        Wave.plots(A1.magsqr().rename('pump'),A2.magsqr().rename('seed'),A3.magsqr().rename('SFG'),
                lines=[0,L/v1,L/v2,L/v3],x='time (ps)',y='power (W)',grid=1,log=1,
                xlim=(L/v1-4*τ,L/v3+4*τ),ylim=(1e-4,1e5),corner='ul',
                save=f'sfg pulse propagation {τ:g}ps {L:g}mm {λ1:g}nm {λ2:g}nm {sell} {η:g}η log')
        Wave.plots(A1.magsqr().rename('pump'),A2.magsqr().rename('seed'),A3.magsqr().rename('SFG'),
                lines=[0,L/v1,L/v2,L/v3],x='time (ps)',y='power (W)',grid=1,log=0,
                xlim=(L/v1-4*τ,L/v3+4*τ),ylim=(-0.1*P1,1.1*P1),corner='ul',
                save=f'sfg pulse propagation {τ:g}ps {L:g}mm {λ1:g}nm {λ2:g}nm {sell} {η:g}η')
    return A1,A2,A3
# @memory.cache
def triplerpulseprop(τ, L, P1, η1=0, η2=0, λ1=1550, sell='ktp', Type='zzz', Λ1=None, Λ2=None, d1=None, d2=None, nres=40, rtol=1e-3, atol=1e-6, plot=False):
    # τ = fwhm pulse duration in ps, L in mm, P = peak power in W, η in %/W/cm²
    λ2,λ3 = λ1/2,λ1/3
    dt,c = τ/nres,0.299792458 # mm/ps = m/ns
    ω1,ω2,ω3 = (1e6*2*np.pi*c/λ for λ in (λ1,λ2,λ3))
    ωmax,ωmin = (1e6*2*np.pi*c/λ for λ in (100,6000)) # min,max λ in nm
    v1,v2,v3 = (c*groupvelocity(λ,sell) for λ in (λ1, λ2, λ3))
    kg1,kg2 = 0 if Λ1 is None else 1e3*2*pi*(1/polingperiod(λ1, λ1, sell, Type) - 1/Λ1), 0 if Λ2 is None else 1e3*2*pi*(1/polingperiod(λ1, λ2, sell, Type) - 1/Λ2)
    g1,g2 = sqrt(0.01*0.01*η1), sqrt(0.01*0.01*η2) # 1/100x for %, 1/100x for cm² to mm²
    d1,d2 = d1 if d1 is not None else lambda z:1, d2 if d2 is not None else lambda z:1
    ts = np.linspace(-100*τ, L/v3+100*τ, 1+int((L/v3+200*τ)//dt))  # ps
    N = len(ts)
    E1 = sqrt(P1 * np.exp(-4*np.log(2)*ts**2/τ**2)) + 0j
    E2,E3 = 0*E1,0*E1
    ωs = 2*np.pi*np.fft.fftfreq(N, dt)
    with np.errstate(invalid='ignore'):
        λ1s,λ2s,λ3s = [1e6*2*np.pi*c / np.clip(ω+ωs,ωmin,ωmax) for ω in (ω1,ω2,ω3)]
        prop1 = -1j*index(λ1s,sell)*(ω1+ωs)/c + 1j*index(λ1,sell)*ω1/c
        prop2 = -1j*index(λ2s,sell)*(ω2+ωs)/c + 1j*index(λ2,sell)*ω2/c
        prop3 = -1j*index(λ3s,sell)*(ω3+ωs)/c + 1j*index(λ3,sell)*ω3/c
    def ode(z, y):
        E1z, E2z, E3z = np.split(y, 3)
        dE1 = np.fft.ifft(prop1*np.fft.fft(E1z))
        dE2 = np.fft.ifft(prop2*np.fft.fft(E2z))
        dE3 = np.fft.ifft(prop3*np.fft.fft(E3z))
        dE1 += 1j*g1 * d1(z/L) * np.conj(E1z)*E2z * exp(1j*kg1*z)       # shg
        dE2 += 1j*g1 * np.conj(d1(z/L)) * E1z**2 * exp(-1j*kg1*z)       # shg
        dE1 += 1j*g2*(λ3/λ1) * d2(z/L) * E2z.conj()*E3z * exp(1j*kg2*z) # sfg
        dE2 += 1j*g2*(λ3/λ2) * d2(z/L) * E1z.conj()*E3z * exp(1j*kg2*z) # sfg
        dE3 += 1j*g2 * np.conj(d2(z/L)) * E1z*E2z * exp(-1j*kg2*z)      # sfg
        return np.concatenate([dE1, dE2, dE3])
    y0 = np.concatenate([E1, E2, E3])
    sol = scipy.integrate.solve_ivp(ode, [0,L], y0, atol=atol, rtol=rtol, t_eval=[0,L], method='DOP853') # RK45, BDF, Radau, LSODA, RK23
    A1,A2,A3 = (Wave(sol.y[i*N:(i+1)*N, -1], ts) for i in range(3))
    A0 = Wave(E1,ts)
    U0,U1,U2,U3 = [A.magsqr().sum()*A.dx() for A in (A0,A1,A2,A3)]
    # print(U0,'pJ U0')
    # print(U1+U2+U3,'pJ U1+U2+U3')
    if plot:
        Wave.plots(A1.magsqr().rename('pump'),A2.magsqr().rename('seed'),A3.magsqr().rename('SFG'),
                lines=[0,L/v1,L/v2,L/v3],x='time (ps)',y='power (W)',grid=1,log=1,
                xlim=(L/v1-4*τ,L/v3+4*τ),ylim=(1e-4,1e5),corner='ul',
                save=f'sfg pulse propagation {τ:g}ps {P1:g}W {L:g}mm {λ1:g}nm {λ2:g}nm {sell} {η1:g}ηshg {η2:g}ηsfg log')
        Wave.plots(A1.magsqr().rename('pump'),A2.magsqr().rename('seed'),A3.magsqr().rename('SFG'),
                lines=[0,L/v1,L/v2,L/v3],x='time (ps)',y='power (W)',grid=1,log=0,
                xlim=(L/v1-4*τ,L/v3+4*τ),ylim=(-0.1*P1,1.1*P1),corner='ul',
                save=f'sfg pulse propagation {τ:g}ps {P1:g}W {L:g}mm {λ1:g}nm {λ2:g}nm {sell} {η1:g}ηshg {η2:g}ηsfg')
    return A0,A1,A2,A3
# @timeit
# @memory.cache
def threewavepulseprop(P1=0, P2=0.01, P3=5000, τ1=50, τ2=50, τ3=0.5, L=1, η=600, λ1=920, λ2=1170, sell='mglnridgewg', Type='zzz', Λ=None, d=None, dt=0.01, rtol=1e-3, atol=1e-6, plot=False):
    # P1,P2,P3 = pulse peak powers in W, τ1,τ2,τ3 = fwhm pulse duration in ps, L in mm, η in %/W/cm²
    c = 0.299792458 # mm/ps = m/ns
    ω1, ω2 = 1e6*2*np.pi*c/λ1, 1e6*2*np.pi*c/λ2
    λ3,ω3 = 1/(1/λ1 + 1/λ2), ω1 + ω2
    v1, v2, v3 = (c*groupvelocity(λ,sell) for λ in (λ1, λ2, λ3))
    kg = 0 if Λ is None else 1e3*2*pi * (1/polingperiod(λ1,λ2,sell,Type) - 1/Λ) # in mm⁻¹
    if Λ:
        assert 0, 'double check that these agree'
        kg = 0 if Λ is None else index(λ3, sell)*ω3/c - index(λ1, sell)*ω1/c - index(λ2, sell)*ω2/c - 2*np.pi*(1e3/Λ) # https://chatgpt.com/c/68acbac0-ae4c-8321-921b-0da64244b171
    g3 = sqrt(0.01*0.01*η) # 1/100x for %, 1/100x for cm² to mm²
    g1,g2 = (λ3/λ1)*g3,(λ3/λ2)*g3
    d = d if d is not None else lambda z:1
    P,τ = max(P1, P2, P3),min(τ1, τ2, τ3)
    ts = wrange(-20*τ, int(L/v3/dt)*dt+20*τ, dt)  # ps
    N = len(ts)
    E1,E2,E3 = [sqrt(P * np.exp(-4*np.log(2)*ts**2/τ**2)) + 0j for P,τ in [(P1,τ1),(P2,τ2),(P3,τ3)]]
    ωs = 2*np.pi*np.fft.fftfreq(N, dt)
    λ1s,λ2s,λ3s = [1e6*2*np.pi*c / (ω + ωs) for ω in (ω1,ω2,ω3)]
    prop1 = -1j*index(λ1s,sell)*(ω1+ωs)/c + 1j*index(λ1,sell)*ω1/c
    prop2 = -1j*index(λ2s,sell)*(ω2+ωs)/c + 1j*index(λ2,sell)*ω2/c
    prop3 = -1j*index(λ3s,sell)*(ω3+ωs)/c + 1j*index(λ3,sell)*ω3/c
    def ode(z, y):
        E1z, E2z, E3z = np.split(y, 3)
        dE1 = np.fft.ifft(prop1*np.fft.fft(E1z))
        dE2 = np.fft.ifft(prop2*np.fft.fft(E2z))
        dE3 = np.fft.ifft(prop3*np.fft.fft(E3z))
        dE1 += 1j*g1 * d(z/L) * E2z.conj()*E3z * np.exp(1j*kg*z)
        dE2 += 1j*g2 * d(z/L) * E1z.conj()*E3z * np.exp(1j*kg*z)
        dE3 += 1j*g3 * conj(d(z/L)) * E1z*E2z * np.exp(-1j*kg*z)
        return np.concatenate([dE1, dE2, dE3])
    y0 = np.concatenate([E1, E2, E3])
    sol = scipy.integrate.solve_ivp(ode, [0,L], y0, atol=atol, rtol=rtol, t_eval=[0,L], method='DOP853') # RK45, BDF, Radau, LSODA, RK23
    # sol = scipy.integrate.solve_ivp(ode, [0,L], y0, atol=atol, rtol=rtol, method='LSODA', t_eval=[0,L], first_step=0.0001, min_step=0.0001, max_step=0.001) # RK45, BDF, Radau, LSODA, RK23
    A1,A2,A3 = (Wave(sol.y[i*N:(i+1)*N, -1], ts) for i in range(3))
    e1,e2,e3 = Wave(E1,ts),Wave(E2,ts),Wave(E3,ts)
    if plot:
        ws = [A1.magsqr().rename('λ1'),A2.magsqr().rename('λ2'),A3.magsqr().rename('λ3')]
        args = dict(lines=[0,L/v1,L/v2,L/v3],x='time (ps)',y='power (W)',grid=1,xlim=(L/v1-4*τ,L/v3+4*τ),corner='ul')
        Wave.plots(*ws,**args,log=1,ylim=(1e-4,1e5),save=f'threewave pulse propagation {τ:g}ps {L:g}mm {λ1:.4g}nm {λ2:.4g}nm {sell} {η:g}η log')
        Wave.plots(*ws,**args,log=0,ylim=(-0.1*P,1.1*P),save=f'threewave pulse propagation {τ:g}ps {L:g}mm {λ1:.4g}nm {λ2:.4g}nm {sell} {η:g}η')
    return A1,A2,A3,e1,e2,e3
def opapulseprop(P3=5e3, P2=0.01, τ3=0.5, τ2=50, L=1, η=600, λ3=515, λ2=1170, sell='mglnridgewg', Type='zzz', Λ=None, d=None, dt=0.01, rtol=1e-3, atol=1e-6, plot=False):
    # OPA/DFG as a special case of three-wave mixing: pump ω3, seed ω2, idler ω1
    λ1,τ1 = 1/(1/λ3 - 1/λ2),max(τ2, τ3)
    A1,A2,A3,e1,e2,e3 = threewavepulseprop(P1=0, P2=P2, P3=P3,τ1=τ1, τ2=τ2, τ3=τ3,L=L, η=η,λ1=λ1, λ2=λ2,sell=sell, Type=Type, Λ=Λ, d=d,dt=dt, rtol=rtol, atol=atol, plot=plot)
    return A1,A2,A3,e2,e3
def dnom(z): return 1
def dhalf(z): return 0.5
def drev(z): return 1-z
def dlin(z): return z
def apod(z): return 1-abs(2*z-1)
def consecd1(z): return z<0.5
def consecd2(z): return z>0.5
def consec2d1(z): return z<0.25
def consec2d2(z): return z>0.25
aa = Wave([0,0.0780203,0.153918,0.226646,0.295616,0.360619,0.421755,0.47936,0.533936,0.586109,0.63662,0.686245,0.7349,0.782245,0.827742,0.870612,0.909812,0.94405,0.971803,0.991343,1],wrange(0,1,0.05))
def dlinmax(z): return aa(z)
def drevmax(z): return aa(1-z)
# @memory.cache
def ηfromd(g,L=2.5):
    aa = Wave([0,0.0780203,0.153918,0.226646,0.295616,0.360619,0.421755,0.47936,0.533936,0.586109,0.63662,0.686245,0.7349,0.782245,0.827742,0.870612,0.909812,0.94405,0.971803,0.991343,1],wrange(0,1,0.05))
    # g = array that defines shg efficiency vs length (sfg efficiency is [1-x for x in g])
    assert all([0<=gi<=1 for gi in g])
    gz = np.linspace(0,1,len(g))
    xs = wrange(0,1,0.01)
    a1s = [aa(Wave(g,gz)(x)) for x in xs]
    a2s = [aa(Wave([1-gi for gi in g],gz)(x)) for x in xs]
    def d1(z): return np.interp(z,xs,a1s)
    def d2(z): return np.interp(z,xs,a2s)
    A0,A1,A2,A3 = triplerpulseprop.func(τ=.22,L=L,P1=1000,η1=27,η2=210,λ1=1556.4,sell='mglnridgewg',d1=d1,d2=d2,nres=40,plot=0)
    η = pulseenergy(A3)/pulseenergy(A0)
    dt = A3.magsqr().fwhm()
    return η,dt

if __name__ == '__main__':
    singlepulseprop(τ=.1,L=5,P=1,λ0=1555,sell='ktp',dt=0.02,plot=1)
    shgpulseprop(τ=.1,L=5,P=1000,η=100,λ0=1555,sell='ktp',dt=0.02,plot=1)
    sfgpulseprop(τ=.2,L=1,P1=1000,P2=1000,η=500,λ1=1550,λ2=775,sell='mglnridgewg',dt=0.02,plot=1)
    triplerpulseprop(τ=.2,L=1,P1=1000,η1=30,η2=250,λ1=1555,sell='mglnridgewg',plot=1)
    threewavepulseprop(plot=1) # opa example by default
    opapulseprop(plot=1)

    import winsound; winsound.Beep(523,500) # C5

    # todo: implement Type
    # todo: In practice, even a few‑percent RMS fluctuation in input peak power can translate into tens of femtoseconds of timing jitter after a centimetre‑length crystal—large enough to degrade pulse‑to‑pulse 
    # coherence in cascaded frequency‑conversion chains or timing‑encoded optical links. Mitigation strategies include operating in the undepleted regime, shortening the crystal so that walk‑off is small 
    # compared with the pulse width, or using aperiodic gratings that flatten Γ(z) over the interaction length so that depletion occurs uniformly across the envelope.
    # less walkoff for type II?
