# ProjectGradys-IME
Arquivos do OMNeT++ para simulação com drones no Aterro do Flamengo. <br />
Colaboradores: <br />
D.Sc. Vítor Gouvêa Andrezo Carneiro (IME) <br />
D.Sc. Bruno José Olivieri de Souza (PUC-Rio) <br />
Thiago Lamenza (PUC-Rio) <br />
Artur Santiago de Oliveira Meneses (IME)
## Introdução
Aqui estão presentes os arquivos e maneiras de uso para realizar a simulação de drone - estação conduzida pela equipe do IME no Aterro do Flamengo, como forma de obter um gráfico de potência recebida na estação por tempo a fim de estudar e modelar efeitos de propagação e canal.
## Instruções para uso
Os arquivos presentes nas pastas src e simulations devem ser carregados em um novo projeto criado no OMNeT++, que seja habilitado em suas configurações para importar e utilizar a biblioteca INET. Ao criar o projeto, os arquivos correspondentes à pasta "src" devem ser carregadas na "src" criada por default. O mesmo procedimento deve ser feito para os arquivos na pasta "simulations".
Os arquivos presentes nas pastas INETDipoleTheoreticalAntenna, INETMobility e INETStaticMobility devem ser usados da seguinte maneira: 
### INETMobility
Dentro da biblioteca INET, o utilizador deve abrir inet4 -> src -> inet -> mobility -> single e copiar os arquivos do Git (INETMobility) para essa pasta, sobreescrevendo os demais que já estavam presentes na biblioteca.
### INETStaticMobility
Dentro da biblioteca INET, o utilizador deve abrir inet4 -> src -> inet -> mobility -> static e copiar os arquivos do Git (INETStaticMobility) para essa pasta, sobreescrevendo os demais que já estavam presentes na biblioteca.
### INETSDipoleTheoreticalAntenna
Dentro da biblioteca INET, o utilizador deve abrir inet4 -> src -> inet -> physicallayer -> antenna e copiar os arquivos do Git (INETDipoleTheoreticalAntenna) para essa pasta, sobreescrevendo os demais que já estavam presentes na biblioteca.
