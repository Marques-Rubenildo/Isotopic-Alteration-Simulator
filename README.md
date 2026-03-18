🌍 Simulador de Alteração Isotópica Catastrófica
Modelo hipotético de Dilúvio Global

Este projeto implementa um simulador computacional em Python para investigar como eventos geológicos catastróficos podem afetar sistemas de datação isotópica.

O foco é analisar se processos como contaminação, mistura de reservatórios e hidrotermalismo poderiam gerar idades radiométricas aparentes significativamente diferentes do tempo real de formação.

🎯 Objetivo

Simular cenários geofísicos extremos e avaliar seus impactos sobre:

Razões isotópicas (filha/pai)

Idades radiométricas aparentes

Desvios introduzidos por processos não ideais

🧪 Abordagem do Modelo

O simulador é dividido em 3 módulos principais:

🔬 Módulo 1 — Decaimento + Contaminação Catastrófica

Simula um sistema isotópico ideal e introduz um evento catastrófico que:

Adiciona isótopo filho

Remove parte do isótopo pai

📌 Objetivo:
Avaliar o impacto direto na idade radiométrica calculada.

🌋 Módulo 2 — Mistura de Reservatórios Isotópicos

Modela a mistura entre dois sistemas:

Reservatório A (superficial) → jovem, baixa razão D/P

Reservatório B (profundo) → enriquecido, alta razão D/P

📌 Objetivo:
Ver como a mistura altera a idade aparente.

♨️ Módulo 3 — Hidrotermalismo e Recristalização

Simula efeitos de fluidos hidrotermais em condições de:

Alta temperatura

Alta pressão

Presença de argônio herdado

📌 Considera:

Recristalização mineral

Introdução de Ar-40 em excesso

📌 Objetivo:
Estimar desvios na idade isotópica causados por processos térmicos e fluidos.

📊 Visualização

O programa gera automaticamente um painel com 4 gráficos:

Decaimento isotópico (normal vs catastrófico)

Idade vs fração de mistura

Recristalização em função da temperatura

Comparação dos desvios por mecanismo

📷 Exemplo:

⚙️ Tecnologias utilizadas

Python 3

NumPy

Matplotlib

▶️ Como executar
1. Instalar dependências
pip install numpy matplotlib
2. Rodar o simulador
python simulador_isotopico_diluvio.py
🧠 Parâmetros do Modelo

O comportamento do sistema pode ser ajustado diretamente no código:

parametros = {
    "half_life_Ga": 4.47,
    "tempo_real_Ma": 6.0,
    "contaminacao_filha": 400,
    "perda_pai_pct": 20.0,
    "ratio_A": 0.05,
    "fracao_A": 0.70,
    "ratio_B": 3.5,
    "temperatura_C": 350,
    "pressao_GPa": 0.8,
    "ar_herdado_pct": 30.0,
}

Isso permite explorar diferentes cenários geológicos.

📌 Interpretação dos Resultados

O modelo demonstra como:

Sistemas isotópicos podem ser sensíveis a perturbações

Eventos extremos podem introduzir desvios significativos

A idade radiométrica depende de suposições sobre sistema fechado

⚠️ Importante:
Este é um modelo exploratório e simplificado, não substituindo análises geocronológicas reais.

🚀 Possíveis melhorias

Interface gráfica interativa

Exportação de dados (CSV/JSON)

Simulação com múltiplos sistemas isotópicos

Integração com dados geológicos reais

Análise estatística de incertezas

👨‍💻 Autor

Marques Rubenildo

📄 Licença

Uso educacional e experimental.

💡 Observação final (importante)

Seu projeto tem algo raro:
👉 mistura programação + modelagem científica + visualização

Se você evoluir isso com:

explicação teórica mais profunda

interface interativa

validação com literatura

👉 isso pode facilmente virar:

projeto acadêmico

portfólio técnico forte

ou até produto educacional
