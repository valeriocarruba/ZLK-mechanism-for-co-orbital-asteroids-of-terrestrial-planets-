import rebound
import os
import numpy as np
import matplotlib.pyplot as plt
from angles import normalize
import concurrent.futures
from astroquery.jplhorizons import Horizons
import astropy.table



# -----------------------------
# Lê os elementos orbitais dos asteroides
# -----------------------------
def ler_elementos_particles(arquivo="Particles.el"):
    with open(arquivo, "r") as f:
        linhas = f.readlines()
    elementos = []
    for linha in linhas:
        partes = linha.strip().split()
        a, e, inc, Omega, omega, M = map(float, partes[:6])
        elementos.append((a, e, np.radians(inc), np.radians(Omega), np.radians(omega), np.radians(M)))
    return elementos

# -----------------------------
# Dados dos planetas
# -----------------------------

bodies = ["199","299","399","301","499","599","699","799","899"]


all_results = astropy.table.vstack([
    Horizons(id=i, location='500@10', epochs=2460800.5).elements()
    for i in bodies
])


# Massas dos planetas e corpos em unidades de massa solar (Msun)
massas = [
    1.6601208254808336e-07,  # Mercúrio
    2.447838287784771e-06,   # Vênus
    3.0034896154502038e-06,  # Terra
    3.694303350091508e-08,  # Moon
    3.2271559174983593e-07,  # Marte
    9.545942479871165e-04,   # Júpiter
    2.8581500081698117e-04,  # Saturno
    4.365793681773934e-05,   # Urano
    5.1503084159155005e-05   # Netuno
]


# -----------------------------
# Tempo de simulação
# -----------------------------
N = 7858
times = np.linspace(0., 0.7857000e5, N)

# -----------------------------
# Simulação de um asteroide
# -----------------------------
def simular_asteroide(args):
    idx, elementos_asteroide = args
    a_ast, e_ast, inc_ast, Omega_ast, omega_ast, M_ast = elementos_asteroide

    sim = rebound.Simulation()
    sim.units = ('AU', 'yr', 'Msun')
    sim.integrator = "bs"

    sim.add(m=0.9999999999950272)  # Sol

    i = 0

    # Planetas internos + Lua
    for row in all_results[0:5]:
        # Get the designation (some might be missing the name in parentheses)
        designation = row['targetname']

        sim.add(m=massas[i], a=row["a"], e=row["e"], inc=np.radians(row["incl"]), Omega=np.radians(row["Omega"]), omega=np.radians(row["w"]), M=np.radians(row["M"]), primary = sim.particles[0] )
        i += 1


    # Asteroide
    sim.add(
        m=0, a=a_ast, e=e_ast, inc=inc_ast, Omega=Omega_ast,
        omega=omega_ast, M=M_ast, primary=sim.particles[0]
    )



    # Planetas externos
    for row in all_results[5:]:
        # Get the designation (some might be missing the name in parentheses)
        designation = row['targetname']

        sim.add(m=massas[i], a=row["a"], e=row["e"], inc=np.radians(row["incl"]), Omega=np.radians(row["Omega"]), omega=np.radians(row["w"]), M=np.radians(row["M"]), primary = sim.particles[0] )
        i += 1


    sim.move_to_com()
    #sim.status()

    e_ast_time = np.zeros(N)
    omega_ast_time = np.zeros(N)
    tv_ast_time = np.zeros(N)
    for i, t in enumerate(times):
        sim.integrate(t)
        orb = sim.particles[6].orbit(primary=sim.particles[0])
        e_ast_time[i] = orb.e
        omega_ast_time[i] = orb.omega*180.0/np.pi
        tv_ast_time[i] = sim.t

    # Plot
    fig2, axs = plt.subplots(2, 1, figsize=(10, 12))

    axs[0].plot(omega_ast_time, e_ast_time, '.')
    axs[0].set_title('Kozai-Resonance')
    axs[0].set_xlabel(r'$\omega$')
    axs[0].set_ylabel(r'$e$ (eccentricity)')

    axs[1].plot(tv_ast_time, omega_ast_time)
    axs[1].set_title('Argument of Periapsis vs Time')
    axs[1].set_xlabel('Time')
    axs[1].set_ylabel(r'$\omega$')

    plt.tight_layout()

    filename = f"figure_{idx:02d}.png"
    plt.savefig(filename)
    plt.close()

    return filename

# -----------------------------
# Execução principal
# -----------------------------
def main():
    print(f"Working directory: {os.getcwd()}")
    asteroides = ler_elementos_particles()
    print(f"{len(asteroides)} asteroides carregados.")

    args_list = list(enumerate(asteroides, start=1))

    from tqdm import tqdm
    with concurrent.futures.ProcessPoolExecutor(max_workers=6) as executor:
        resultados = list(tqdm(executor.map(simular_asteroide, args_list), total=len(args_list)))

    print("Simulações concluídas. Arquivos gerados:")
    for filename in resultados:
        print(filename)

if __name__ == "__main__":
    main()
