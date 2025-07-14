# Cálculo de Splines Cúbicos em C  
**Projeto 2 - Métodos Numéricos e Algoritmos de Programação**

Este repositório documenta o desenvolvimento do Projeto 2, realizado como atividade avaliativa para as disciplinas de *Métodos Numéricos* e *Algoritmos de Programação*. O trabalho consiste na implementação e análise de dois métodos de interpolação: **Spline Cúbico Natural** e **Spline Cúbico Fixado**, utilizando a linguagem C.

O projeto integra conhecimentos das duas áreas:

- **Do ponto de vista de Algoritmos**, envolve a manipulação de arquivos de entrada (`entrada.txt`) e saída (`saida.txt`), a estruturação de dados em arrays para armazenar os pontos e coeficientes, e a modularização do código em funções coesas e reutilizáveis.

- **Sob a ótica de Métodos Numéricos**, o foco está na aplicação da teoria de interpolação, na resolução do sistema de equações tridiagonal para encontrar os coeficientes do spline, e na análise da precisão do método através do cálculo de integrais, derivadas e seus respectivos erros em comparação com uma função de referência.

O objetivo final é demonstrar a aplicação prática de conceitos de programação e cálculo numérico na solução de um problema clássico de interpolação de dados.

---

## ✨ Funcionalidades

- **Dois Tipos de Spline**  
  Implementação tanto do Spline Natural quanto do Spline Fixado.

- **Entrada de Dados via Arquivo**  
  Lê o número de pontos e suas coordenadas (x, y) do arquivo `entrada.txt`.

- **Saída Detalhada**  
  Salva todos os resultados, de forma organizada, no arquivo `saida.txt`.

- **Visualização dos Polinômios**  
  Exibe a equação do polinômio cúbico para cada intervalo entre os pontos.

- **Análise Numérica**
  - Calcula a integral definida do spline em todo o intervalo.
  - Calcula a primeira e a segunda derivada do spline em um ponto pré-definido (`x = 0.5`).
  - Compara os resultados com os valores reais da função `cos(π * x)` e exibe o erro absoluto.

---

## 📂 Estrutura do Projeto

```
.
├── spline_fixado.c       # Código-fonte para o Spline Cúbico Fixado
├── spline_natural.c      # Código-fonte para o Spline Cúbico Natural
├── entrada.txt           # Arquivo de exemplo para entrada de dados
└── README.md             # Este arquivo
```

---

## 🚀 Como Usar

Siga os passos abaixo para compilar e executar os programas.

### 🔧 Pré-requisitos

- Um compilador C (como o GCC).

### 1. Preparar o Arquivo de Entrada

Crie um arquivo chamado `entrada.txt` no mesmo diretório dos códigos.  
O arquivo deve ter o seguinte formato:

- A primeira linha deve conter um único número inteiro, que é a quantidade de pontos `N`.
- As `N` linhas seguintes devem conter os pares de coordenadas `x` e `y`, separados por um espaço.

**Exemplo de `entrada.txt` com 5 pontos:**

```
5
0.0 1.0000000000
0.25 0.7071067812
0.5 0.0000000000
0.75 -0.7071067812
1.0 -1.0000000000
```

### 2. Compilar os Códigos

Abra um terminal e use os seguintes comandos para compilar cada programa.

> **Nota:** A flag `-lm` é necessária para vincular a biblioteca matemática (`math.h`).

**Para compilar o Spline Fixado:**

```bash
gcc spline_fixado.c -o spline_fixado -lm
```

**Para compilar o Spline Natural:**

```bash
gcc spline_natural.c -o spline_natural -lm
```

### 3. Executar o Programa

Após a compilação, execute o programa desejado.  
Certifique-se de que o arquivo `entrada.txt` está presente.

**Para executar o Spline Fixado:**

```bash
./spline_fixado
```

**Para executar o Spline Natural:**

```bash
./spline_natural
```

### 4. Verificar a Saída

Após a execução, uma mensagem de conclusão aparecerá no terminal.  
Um arquivo chamado `saida.txt` será criado (ou sobrescrito) com os resultados detalhados, incluindo:

- Um cabeçalho confirmando o tipo de spline e os dados lidos.
- A lista completa das equações polinomiais para cada intervalo.
- A análise final com os valores da integral e das derivadas, junto com os respectivos erros.

---
