Pode utilizar o seguinte comando bash aliado ao python pra verificar se a energia convergiu em um diretorio de outputs do QE.

```bash
!cd ../scf_025_NOSYM.out/; for file in *; do \
    awk '/^!.*total energy/ { print FILENAME ": " $5; found=1 } \
         END { if (!found) print FILENAME ": não convergiu" }' "$file"; \
done
```
Atribuimos o resultado acima para `dados_energia`  

Em seguida criamos uma lista: `dados_energia = [ valor.split(': ',) for valor in dados_energia]` para entrão criar um DataFrame.


df = pd.DataFrame(dados_energia, columns=["Arquivo", "Energia Total"])
df['Energia Total'] = df['Energia Total'].astype(float)

> Supondo que seu DataFrame se chame df
df[['a', 'ratio_ca', 'num_estrutura']] = df['Arquivo'].str.extract(
    r'scf-([0-9.]+)-([0-9.]+)-[0-9.]+-([0-9]+)\.out'
)

> Convertendo para os tipos numéricos apropriados
df['a'] = df['a'].astype(float)
df['ratio_ca'] = df['ratio_ca'].astype(float)
df['num_estrutura'] = df['num_estrutura'].astype(int)


> Contas em que a energia não convergiu.
df[
    df['Energia Total'].isna()
]

> Agrupando por estrutura
df1 = df[df['num_estrutura'] == 1]
df2 = df[df['num_estrutura'] == 2]
df3 = df[df['num_estrutura'] == 3]
df4 = df[df['num_estrutura'] == 4]
df5 = df[df['num_estrutura'] == 5]

lista_dfs = [df1, df2, df3,df4,df5] 
