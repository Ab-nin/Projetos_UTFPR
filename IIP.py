import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

# Banco de dados com nome de moléculas 
mol_db = [
    '1-Butene', 'Acetone', 'Air', 'Ammonia', 'Argon', 'Benzene', 'CarbonDioxide',
    'CarbonMonoxide', 'CarbonylSulfide', 'CycloHexane', 'CycloPropane', 'Cyclopentane',
    'D4', 'D5', 'D6', 'Deuterium', 'Dichloroethane', 'DiethylEther', 'DimethylCarbonate',
    'DimethylEther', 'Ethane', 'Ethanol', 'EthylBenzene', 'Ethylene', 'EthyleneOxide',
    'Fluorine', 'HFE143m', 'HeavyWater', 'Helium', 'Hydrogen', 'HydrogenChloride',
    'HydrogenSulfide', 'IsoButane', 'IsoButene', 'Isohexane', 'Isopentane', 'Krypton',
    'MD2M', 'MD3M', 'MD4M', 'MDM', 'MM', 'Methane', 'Methanol', 'MethylLinoleate',
    'MethylLinolenate', 'MethylOleate', 'MethylPalmitate', 'MethylStearate', 'Neon',
    'Neopentane', 'Nitrogen', 'NitrousOxide', 'Novec649', 'OrthoDeuterium', 'OrthoHydrogen',
    'Oxygen', 'ParaDeuterium', 'ParaHydrogen', 'Propylene', 'Propyne', 'R11', 'R113',
    'R114', 'R115', 'R116', 'R12', 'R123', 'R1233zd(E)', 'R1234yf', 'R1234ze(E)',
    'R1234ze(Z)', 'R124', 'R1243zf', 'R125', 'R13', 'R1336mzz(E)', 'R134a', 'R13I1',
    'R14', 'R141b', 'R142b', 'R143a', 'R152A', 'R161', 'R21', 'R218', 'R22', 'R227EA',
    'R23', 'R236EA', 'R236FA', 'R245ca', 'R245fa', 'R32', 'R365MFC', 'R40', 'R404A',
    'R407C', 'R41', 'R410A', 'R507A', 'RC318', 'SES36', 'SulfurDioxide', 'SulfurHexafluoride',
    'Toluene', 'Water', 'Xenon', 'cis-2-Butene'
]

def plot_poco(modelo, nome_mol_a, sigma_a, nome_mol_b, sigma_b, lambb, k):
    if modelo == 'Poço Quadrado':
        mol_a = nome_mol_a
        T_ca = PropsSI('Tcrit', mol_a)
        ep_s_k_a = T_ca / 1.25
        mol_b = nome_mol_b
        T_cb = PropsSI('Tcrit', mol_b)
        ep_s_k_b = T_cb / 1.25
        distt = (sigma_a + sigma_b) / 2
        comp_poco = distt * lambb
        epi_12 = (((ep_s_k_a * ep_s_k_b) ** (1/2)) * (1 - k))
        
        # Mostrar epsilon_12
        st.latex(r'\epsilon_{12} = %.2f \, K' % epi_12)

        # Coordenadas da plotagem
        x_coords = [distt, distt, comp_poco, comp_poco, 2.0]
        y_coords = [200, -epi_12, -epi_12, 0, 0]

        fig, ax = plt.subplots(dpi=100)
        ax.plot(x_coords, y_coords, '-', linewidth=1.5)
        ax.set_xlim(0.3, 2.0)
        ax.set_ylim(-500, 200)
        ax.set_xlabel('Distância r (nm)')
        ax.set_ylabel('Energia u/k (K)')
        ax.set_title('Potencial de Poço Quadrado')
        ax.grid(True)
        
        return fig

    elif modelo == 'Lennard-Jones':
        # Encontrando parâmetros
        sigma_ab = (sigma_a + sigma_b) / 2
        r_min = (2 ** (1/6)) * sigma_ab

        # Descobrindo os epsilons
        mol_a = nome_mol_a
        T_ca = PropsSI('Tcrit', mol_a)
        ep_s_k_a = T_ca / 1.25
        mol_b = nome_mol_b
        T_cb = PropsSI('Tcrit', mol_b)
        ep_s_k_b = T_cb / 1.25
        epi_12 = (((ep_s_k_a * ep_s_k_b) ** (1/2)) * (1 - k))
        
        # Mostrar parâmetros
        st.latex(r'\epsilon_{12} = %.2f \, K' % epi_12)
        st.latex(r'r_{min} = %.3f \, nm' % r_min)

        # Criar array para os valores de r
        r_valores = np.linspace(0.8*sigma_ab, 3*sigma_ab, 2000)

        # Calcular potencial LJ (u/k em K)
        Vr_valores = 4 * epi_12 * ((sigma_ab / r_valores)**12 - (sigma_ab / r_valores)**6)

        # Plotar
        fig, ax = plt.subplots(figsize = (6,4))
        ax.plot(r_valores, Vr_valores, '-', linewidth=1.5)
        ax.set_xlabel("Distância r (nm)")
        ax.set_ylabel("Energia potencial u/k (K)")
        ax.set_title("Potencial de Lennard-Jones")
        ax.grid(True)
        ax.set_ylim(-2*epi_12, 100)
        
        return fig

# Configuração da página
st.set_page_config(page_title="Potenciais Moleculares", layout="wide")

# Título
st.title("Gráfico de Potenciais Intermoleculares")
st.markdown(
    """
        Este aplicativo permite visualizar os potenciais de interação molecular entre pares de moléculas 
        utilizandoos modelos de Poço Quadrado (Square Well) e Lennard-Jones.

        ### Modo de usar:

        1. Selecione o modelo desejado na barra lateral;
        2. Escolha as moléculas A e B do banco de dados;
        3. Ajuste os parâmetros usando os controles deslizantes;

    """
)

# Sidebar com controles
st.sidebar.header("Parâmetros")

modelo = st.sidebar.selectbox(
    "Modelo:",
    ['Poço Quadrado', 'Lennard-Jones']
)

st.sidebar.subheader("Molécula A")
nome_mol_a = st.sidebar.selectbox(
    "Molécula A:",
    mol_db,
    index=mol_db.index('Methane')
)

sigma_a = st.sidebar.slider(
    "σ_a (nm)",
    min_value=0.2,
    max_value=0.6,
    value=0.36,
    step=0.01
)

st.sidebar.subheader("Molécula B")
nome_mol_b = st.sidebar.selectbox(
    "Molécula B:",
    mol_db,
    index=mol_db.index('Benzene')
)

sigma_b = st.sidebar.slider(
    "σ_b (nm)",
    min_value=0.2,
    max_value=0.6,
    value=0.52,
    step=0.01
)

st.sidebar.subheader("Parâmetros de Interação")
lambb = st.sidebar.slider(
    "Potencial do Poço-Quadrado, λ",
    min_value=0.5,
    max_value=2.5,
    value=1.5,
    step=0.05
)

k = st.sidebar.slider(
    "Parâmetro de interação binário, k",
    min_value=-2.5,
    max_value=3.5,
    value=0.2,
    step=0.05
)

# Gerar e mostrar o gráfico
try:
    fig = plot_poco(modelo, nome_mol_a, sigma_a, nome_mol_b, sigma_b, lambb, k)
    st.pyplot(fig)
except Exception as e:
    st.error(f"Erro ao gerar o gráfico: {str(e)}")

# Informações adicionais
with st.expander("ℹ️ Sobre os Modelos"):
    st.markdown("""
    **Poço Quadrado**: Modelo simplificado que representa a interação molecular como um poço de potencial retangular.
    
    **Lennard-Jones**: Modelo mais realista que descreve a interação entre pares de moléculas neutras, 
    combinando repulsão de curto alcance e atração de longo alcance.
    """)