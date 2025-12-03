#!/usr/bin/env python
import pandas as pd
import csv
import matplotlib.pyplot as plt
import seaborn as sns


class Enlace:
    def __init__(self, acceptor, donor, frac):
        self.acceptor = acceptor
        self.donor = donor
        self.frac = float(frac)  # Convert to float


enlaces = []
data1 = '6UML-POM.out'   #It comes from analysis_cpptraj.sh            


def renombrar_residuo(residuo): 
    resname = residuo.split('@',1)[0]
    return resname.replace('_',' ')
    
    
def intercambio(enlace):
    acceptor_number = int(enlace.acceptor.split(' ')[1])
    donor_number = int(enlace.donor.split (' ')[1])
    aux_acceptor = enlace.acceptor
    aux_donor = enlace.donor
    if acceptor_number not in range (1, 376):
      enlace.acceptor = aux_donor
      enlace.donor = aux_acceptor
    enlace.acceptor = suma_47(enlace.acceptor)
    enlace.donor = suma_26(enlace.donor)  
    return enlace

''' '''      

def suma_47(residuo):
    acceptor_number = int(residuo.split(' ')[1])
    return residuo.split(' ')[0] + " " + str(acceptor_number + 47)


def suma_26(residuo):
    donor_number = int(residuo.split (' ')[1])
    return residuo.split(' ')[0] + " " + str(donor_number + 26)


def agrupar_y_sumar(frac_list):
    agrupados = {}
    for enlace in frac_list:
        llave = (enlace.acceptor, enlace.donor)
        if llave not in agrupados:
            agrupados[llave] = enlace.frac
        else:
            agrupados[llave] += enlace.frac
    
    resultado = [Enlace(acc, don, frac) for (acc, don), frac in agrupados.items()]
    return resultado

def normalizar_frac(frac_list):
    max_frac = max(enlace.frac for enlace in frac_list)
    min_frac = min(enlace.frac for enlace in frac_list)
    for enlace in frac_list:
      enlace.frac = (enlace.frac - min_frac) / (max_frac - min_frac)
    return frac_list



with open(data1, "r") as f:
    next(f)
    for line in f:
        parts = line.split()
        if len(parts) >= 5:
            enlace = Enlace(renombrar_residuo(parts[0]), renombrar_residuo(parts[2]), parts[4])
            enlace = intercambio(enlace)
            enlaces.append(enlace)

enlaces = agrupar_y_sumar(enlaces)
enlaces = normalizar_frac(enlaces)

# Guardar en CSV
with open("enlaces.csv", "w", newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["acceptor", "donor", "frac"])  
    for enlace in enlaces:
        writer.writerow([enlace.acceptor, enlace.donor, enlace.frac])
      
        
df = pd.DataFrame([(e.acceptor, e.donor, e.frac) for e in enlaces], columns=['acceptor', 'donor', 'frac'])
pivot_df = df.pivot(index='acceptor', columns='donor', values='frac')

def ordenar_por_residuo(labels):
    # Ordena las etiquetas extrayendo el numero, asumiendo formato "XXX NNN"
    return sorted(labels, key=lambda x: int(x.split()[1]))

# Ordenar indices y columnas del DataFrame pivote
pivot_df = pivot_df.reindex(ordenar_por_residuo(pivot_df.index), axis=0)
pivot_df = pivot_df.reindex(ordenar_por_residuo(pivot_df.columns), axis=1)

plt.figure(figsize=(12, 12))
ax = sns.heatmap(pivot_df, cmap='Blues',linewidths=0.1, linecolor='gray', xticklabels=True, yticklabels=True) 
plt.title('H-bonds',  fontsize=22)
plt.xlabel('SALL4',  fontsize=18)
plt.ylabel('CRBN',  fontsize=18)
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)

new_x_labels = [label.replace(' ', '\n') for label in pivot_df.columns]
ax.set_xticks([x + 0.5 for x in range(len(new_x_labels))])
ax.set_xticklabels(new_x_labels, rotation=0, ha='center')

plt.tight_layout()
plt.show()
plt.savefig('Hbonds-6UML-POM.png')
