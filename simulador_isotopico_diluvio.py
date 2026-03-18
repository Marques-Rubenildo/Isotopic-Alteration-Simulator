"""
Simulador de Alteração Isotópica Catastrófica
Modelo hipotético de Dilúvio Global

Testa se um evento geológico catastrófico poderia produzir
proporções isotópicas equivalentes às observadas em rochas "antigas".

Requer: numpy, matplotlib
Instalar: pip install numpy matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ─── Constantes ────────────────────────────────────────────────
LN2 = np.log(2)

# ─── Funções base ──────────────────────────────────────────────

def lambda_decay(half_life_Ga):
    """Constante de decaimento a partir da meia-vida em bilhões de anos."""
    return LN2 / (half_life_Ga * 1e9)

def calc_age(ratio_D_P, lam):
    """Idade radiométrica aparente a partir da razão filha/pai."""
    if ratio_D_P <= 0:
        return 0.0
    return np.log(1 + ratio_D_P) / lam

def format_age(years):
    """Formata anos em Ga / Ma legíveis."""
    if years >= 1e9:
        return f"{years/1e9:.3f} Ga"
    elif years >= 1e6:
        return f"{years/1e6:.1f} Ma"
    else:
        return f"{years:.0f} a"


# ══════════════════════════════════════════════════════════════════
# MÓDULO 1 — Decaimento padrão + contaminação catastrófica
# ══════════════════════════════════════════════════════════════════

def modulo_decaimento(
    parent_inicial=1000,
    half_life_Ga=4.5,
    tempo_real_Ma=6.0,
    contaminacao_filha=400,
    perda_pai_pct=20.0,
    tempo_max_Ga=5.0,
    n_pontos=500,
):
    """
    Simula o efeito de uma contaminação catastrófica sobre a razão isotópica.

    Parâmetros
    ----------
    parent_inicial     : quantidade inicial do isótopo pai
    half_life_Ga       : meia-vida em bilhões de anos
    tempo_real_Ma      : tempo real decorrido desde a formação (Ma)
    contaminacao_filha : isótopo filho adicionado pelo evento catastrófico
    perda_pai_pct      : percentual do isótopo pai removido pelo evento
    """
    lam = lambda_decay(half_life_Ga)
    tempo_real = tempo_real_Ma * 1e6

    t = np.linspace(0, tempo_max_Ga * 1e9, n_pontos)
    parent_normal = parent_inicial * np.exp(-lam * t)
    daughter_normal = parent_inicial - parent_normal

    # Aplica contaminação a partir do evento
    evento_idx = np.searchsorted(t, tempo_real)
    daughter_cat = daughter_normal.copy()
    parent_cat = parent_normal.copy()
    daughter_cat[evento_idx:] += contaminacao_filha
    parent_cat[evento_idx:] *= (1 - perda_pai_pct / 100)

    # Cálculo das idades aparentes no momento do evento
    P_real = parent_inicial * np.exp(-lam * tempo_real)
    D_real = parent_inicial - P_real
    ratio_normal = D_real / P_real if P_real > 0 else 0

    D_cat = D_real + contaminacao_filha
    P_cat = P_real * (1 - perda_pai_pct / 100)
    ratio_cat = D_cat / P_cat if P_cat > 0 else 0

    age_normal = calc_age(ratio_normal, lam)
    age_cat = calc_age(ratio_cat, lam)

    print("\n═══ MÓDULO 1 — Decaimento + Contaminação Catastrófica ═══")
    print(f"  Tempo real do evento   : {format_age(tempo_real)}")
    print(f"  Idade aparente padrão  : {format_age(age_normal)}")
    print(f"  Idade aparente catástrofe: {format_age(age_cat)}")
    print(f"  Desvio introduzido     : +{format_age(age_cat - age_normal)}")

    return t, parent_normal, daughter_normal, daughter_cat, parent_cat, tempo_real


# ══════════════════════════════════════════════════════════════════
# MÓDULO 2 — Mistura de reservatórios isotópicos
# ══════════════════════════════════════════════════════════════════

def modulo_mistura(
    ratio_A=0.05,
    fracao_A=0.70,
    ratio_B=3.5,
    half_life_Ga=4.5,
):
    """
    Simula a mistura entre dois reservatórios isotópicos distintos.
    Reservatório A = crosta superficial (jovem, razão baixa)
    Reservatório B = manto profundo (sistema antigo ou enriquecido)
    """
    fracao_B = 1 - fracao_A
    ratio_mix = ratio_A * fracao_A + ratio_B * fracao_B

    lam = lambda_dev = lambda_decay(half_life_Ga)
    age_A = calc_age(ratio_A, lam)
    age_B = calc_age(ratio_B, lam)
    age_mix = calc_age(ratio_mix, lam)

    fracoes = np.linspace(0, 1, 200)
    ages_vs_fracao = [
        calc_age(ratio_A * f + ratio_B * (1 - f), lam)
        for f in fracoes
    ]

    print("\n═══ MÓDULO 2 — Mistura de Reservatórios ═══")
    print(f"  Razão D/P reservatório A  : {ratio_A:.3f}  → {format_age(age_A)}")
    print(f"  Razão D/P reservatório B  : {ratio_B:.3f}  → {format_age(age_B)}")
    print(f"  Fração A na mistura       : {fracao_A*100:.0f}%")
    print(f"  Razão D/P da mistura      : {ratio_mix:.3f}")
    print(f"  Idade aparente da mistura : {format_age(age_mix)}")

    return fracoes, ages_vs_fracao, ratio_mix, age_mix


# ══════════════════════════════════════════════════════════════════
# MÓDULO 3 — Hidrotermalismo e recristalização
# ══════════════════════════════════════════════════════════════════

def modulo_hidrotermal(
    temperatura_C=350,
    duracao_anos=5000,
    pressao_GPa=0.8,
    ar_herdado_pct=30.0,
    half_life_Ga=4.5,
):
    """
    Estima o desvio na idade radiométrica introduzido por
    circulação hidrotermal catastrófica (modelo K-Ar / Ar-Ar).

    O excesso de Argônio herdado é um mecanismo real documentado
    na literatura geológica que pode produzir idades aparentes
    muito superiores à idade real de formação.
    """
    # Fator de recristalização (0–1) dependente de T e P
    fator_temp = min(1.0, (temperatura_C - 100) / 600)
    fator_press = min(1.0, pressao_GPa / 3.0)
    recristalizacao = fator_temp * 0.70 + fator_press * 0.30

    # Excesso de Ar-40 introduzido pelo evento
    # (proporcional à fração herdada e à duração)
    ar_base = ar_herdado_pct / 100
    ar_shift = ar_base * 800e6       # até ~800 Ma de desvio por Ar herdado
    recryst_shift = recristalizacao * 600e6  # até ~600 Ma por recristalização
    total_shift = 50e6 + ar_shift + recryst_shift   # base mínima de 50 Ma

    temps = np.arange(100, 750, 50)
    recryst_por_temp = [min(100, ((t - 100) / 600) * 70 + fator_press * 30) for t in temps]

    print("\n═══ MÓDULO 3 — Hidrotermalismo e Recristalização ═══")
    print(f"  Temperatura do fluido     : {temperatura_C}°C")
    print(f"  Pressão                   : {pressao_GPa} GPa")
    print(f"  Duração do evento         : {duracao_anos:,} anos")
    print(f"  Ar herdado (excesso)      : {ar_herdado_pct:.0f}%")
    print(f"  Grau de recristalização   : {recristalizacao*100:.1f}%")
    print(f"  Desvio total na idade     : +{format_age(total_shift)}")

    return temps, recryst_por_temp, total_shift


# ══════════════════════════════════════════════════════════════════
# VISUALIZAÇÃO
# ══════════════════════════════════════════════════════════════════

def plotar_resultados(params=None):
    if params is None:
        params = {}

    p = {
        "parent_inicial": 1000,
        "half_life_Ga": 4.5,
        "tempo_real_Ma": 6.0,
        "contaminacao_filha": 400,
        "perda_pai_pct": 20.0,
        "ratio_A": 0.05,
        "fracao_A": 0.70,
        "ratio_B": 3.5,
        "temperatura_C": 350,
        "duracao_anos": 5000,
        "pressao_GPa": 0.8,
        "ar_herdado_pct": 30.0,
    }
    p.update(params)

    # ─── Executar módulos ───────────────────────────────────────
    t, pn, dn, dc, pc, t_real = modulo_decaimento(
        parent_inicial=p["parent_inicial"],
        half_life_Ga=p["half_life_Ga"],
        tempo_real_Ma=p["tempo_real_Ma"],
        contaminacao_filha=p["contaminacao_filha"],
        perda_pai_pct=p["perda_pai_pct"],
    )

    fracoes, ages_mix, ratio_mix, age_mix = modulo_mistura(
        ratio_A=p["ratio_A"],
        fracao_A=p["fracao_A"],
        ratio_B=p["ratio_B"],
        half_life_Ga=p["half_life_Ga"],
    )

    temps, recryst, total_shift = modulo_hidrotermal(
        temperatura_C=p["temperatura_C"],
        duracao_anos=p["duracao_anos"],
        pressao_GPa=p["pressao_GPa"],
        ar_herdado_pct=p["ar_herdado_pct"],
        half_life_Ga=p["half_life_Ga"],
    )

    # ─── Figura ────────────────────────────────────────────────
    fig = plt.figure(figsize=(16, 12), facecolor="#0f0f14")
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.35)

    cor_fundo = "#0f0f14"
    cor_painel = "#1a1a24"
    cor_texto = "#e8e8f0"
    cor_muted = "#888899"
    CORES = ["#4a9eff", "#4cdf87", "#ff6b4a", "#f0c040"]

    def estilizar_eixo(ax, titulo):
        ax.set_facecolor(cor_painel)
        ax.tick_params(colors=cor_muted, labelsize=9)
        for spine in ax.spines.values():
            spine.set_edgecolor("#333344")
        ax.set_title(titulo, color=cor_texto, fontsize=11, pad=10)
        ax.xaxis.label.set_color(cor_muted)
        ax.yaxis.label.set_color(cor_muted)

    # ── Gráfico 1: Decaimento + contaminação ───────────────────
    ax1 = fig.add_subplot(gs[0, 0])
    t_Ga = t / 1e9
    ax1.plot(t_Ga, pn, color=CORES[0], linewidth=1.5, label="Pai (padrão)")
    ax1.plot(t_Ga, dn, color=CORES[1], linewidth=1.5, label="Filho (padrão)")
    ax1.plot(t_Ga, dc, color=CORES[2], linewidth=1.5, linestyle="--", label="Filho pós-catástrofe")
    ax1.plot(t_Ga, pc, color=CORES[3], linewidth=1.5, linestyle="--", label="Pai pós-catástrofe")
    ax1.axvline(t_real / 1e9, color="#ffffff", linewidth=0.8, linestyle=":", alpha=0.5)
    ax1.set_xlabel("Tempo (Ga)")
    ax1.set_ylabel("Quantidade (u.a.)")
    ax1.legend(fontsize=8, facecolor=cor_painel, edgecolor="#333344", labelcolor=cor_texto)
    estilizar_eixo(ax1, "Módulo 1 — Decaimento + contaminação")

    # ── Gráfico 2: Razão D/P vs fração da mistura ──────────────
    ax2 = fig.add_subplot(gs[0, 1])
    ages_Ga = [a / 1e9 for a in ages_mix]
    ax2.plot(fracoes * 100, ages_Ga, color=CORES[0], linewidth=2)
    ax2.axhline(age_mix / 1e9, color=CORES[2], linewidth=1, linestyle="--", label=f"Mistura atual: {format_age(age_mix)}")
    ax2.set_xlabel("% do reservatório A na mistura")
    ax2.set_ylabel("Idade aparente (Ga)")
    ax2.legend(fontsize=8, facecolor=cor_painel, edgecolor="#333344", labelcolor=cor_texto)
    estilizar_eixo(ax2, "Módulo 2 — Mistura de reservatórios")

    # ── Gráfico 3: Recristalização por temperatura ──────────────
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.bar(temps, recryst, width=40, color=CORES[0], alpha=0.85)
    ax3.axhline(p["ar_herdado_pct"], color=CORES[2], linewidth=1.5, linestyle="--",
                label=f"Ar herdado: {p['ar_herdado_pct']:.0f}%")
    ax3.set_xlabel("Temperatura (°C)")
    ax3.set_ylabel("% recristalização estimada")
    ax3.set_ylim(0, 110)
    ax3.legend(fontsize=8, facecolor=cor_painel, edgecolor="#333344", labelcolor=cor_texto)
    estilizar_eixo(ax3, "Módulo 3 — Recristalização hidrotermal")

    # ── Gráfico 4: Resumo dos desvios ──────────────────────────
    ax4 = fig.add_subplot(gs[1, 1])
    mecanismos = ["Contaminação\nfilha", "Perda\ndo pai", "Mistura de\nreservatórios", "Ar herdado\n(hidrotermal)"]
    lam = lambda_decay(p["half_life_Ga"])
    P0 = p["parent_inicial"]
    t_r = p["tempo_real_Ma"] * 1e6
    P_r = P0 * np.exp(-lam * t_r)
    D_r = P0 - P_r
    age_base = calc_age(D_r / P_r if P_r > 0 else 0, lam)

    D1 = D_r + p["contaminacao_filha"]
    shift1 = (calc_age(D1 / P_r, lam) - age_base) / 1e6 if P_r > 0 else 0

    P2 = P_r * (1 - p["perda_pai_pct"] / 100)
    shift2 = (calc_age(D_r / P2, lam) - age_base) / 1e6 if P2 > 0 else 0

    shift3 = (age_mix - calc_age(p["ratio_A"], lam)) / 1e6
    shift4 = p["ar_herdado_pct"] / 100 * 800

    desvios = [max(0, shift1), max(0, shift2), max(0, shift3), max(0, shift4)]
    bars = ax4.bar(mecanismos, desvios, color=[CORES[0], CORES[1], CORES[2], CORES[3]], alpha=0.85)
    ax4.set_ylabel("Desvio na idade (Ma)")
    for bar, val in zip(bars, desvios):
        ax4.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 5,
                 f"{val:.0f} Ma", ha="center", va="bottom", fontsize=8, color=cor_texto)
    estilizar_eixo(ax4, "Resumo — Desvio por mecanismo (Ma)")

    # ── Título principal ────────────────────────────────────────
    fig.suptitle(
        "Simulador de Alteração Isotópica Catastrófica\nModelo hipotético de Dilúvio Global",
        color=cor_texto, fontsize=14, y=0.98
    )

    plt.savefig("simulador_diluvio.png", dpi=150, bbox_inches="tight", facecolor=cor_fundo)
    print("\n  Gráfico salvo em: simulador_diluvio.png")
    plt.show()


# ══════════════════════════════════════════════════════════════════
# PONTO DE ENTRADA
# ══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("=" * 60)
    print("  SIMULADOR ISOTÓPICO — MODELO DE DILÚVIO GLOBAL")
    print("=" * 60)

    # ── Edite estes parâmetros para explorar o modelo ──────────
    parametros = {
        # Sistema isotópico (ex.: U-Pb ≈ 4.47 Ga; K-Ar ≈ 1.25 Ga; Rb-Sr ≈ 48.8 Ga)
        "parent_inicial": 1000,
        "half_life_Ga": 4.47,       # meia-vida em bilhões de anos
        "tempo_real_Ma": 6.0,       # tempo REAL do evento (Ma)

        # Contaminação catastrófica
        "contaminacao_filha": 400,  # isótopo filho injetado
        "perda_pai_pct": 20.0,      # % do isótopo pai removido

        # Mistura de reservatórios
        "ratio_A": 0.05,            # razão D/P crosta superficial
        "fracao_A": 0.70,           # fração da crosta na mistura
        "ratio_B": 3.5,             # razão D/P manto profundo

        # Hidrotermalismo
        "temperatura_C": 350,       # temperatura do fluido
        "duracao_anos": 5000,       # duração do evento
        "pressao_GPa": 0.8,         # pressão em GPa
        "ar_herdado_pct": 30.0,     # % de argônio herdado (excesso)
    }

    plotar_resultados(parametros)
