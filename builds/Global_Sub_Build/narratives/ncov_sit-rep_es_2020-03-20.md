---
title: Análisis genómico de la propagación del COVID-19. Informe de la situación hasta el 2020-03-20.
authors:
  - Emma Hodcroft
  - Nicola Müller
  - Cassia Wagner
  - Misja Ilcisin
  - James Hadfield
  - Sidney M. Bell
  - Richard Neher
  - Trevor Bedford
authorLinks:
  - https://neherlab.org/emma-hodcroft.html
  - https://bedford.io/team/nicola-mueller/
  - https://bedford.io/team/cassia-wagner/
  - https://bedford.io/team/misja-ilcisin/
  - https://bedford.io/team/james-hadfield/
  - https://twitter.com/sidneymbell
  - https://neherlab.org/richard-neher.html
  - https://bedford.io/team/trevor-bedford/
affiliations: "Fred Hutch, Seattle, USA; Biozentrum, Basel, Switzerland; CZI, CA, USA"
translators:
  - Roy Costilla
  - Miguel I. Paredes
translatorLinks:
  - https://researchers.uq.edu.au/researcher/18392
  - https://twitter.com/miguelp1120
date: "2020 March 20"
dataset: "https://nextstrain.org/ncov/2020-03-20?legend=closed&d=map&legend=closed"
abstract: "Este reporte utiliza datos genómicos públicos para el seguimiento de la propagación del COVID-19. Los reportes son actualizados semanalmente."
---
<!-- Translators: Only text after : in the above ^ needs to be translated -->
<!-- Comment tags like these do not need to be translated, they are only to help you! -->
<!-- Ensure that links always end in a 'letter' (. counts) If some kind of text doesn't follow them, it breaks the slide. -->
<!-- numbers can be tagged ilke this: 161</tag> - this is just for us to help find them to update! Just leave in the </tag> bit. -->
<!-- Acentos: 	áéíóú -->

<!-- This is left-side text 1 -->
# [Contenidos](https://nextstrain.org/ncov/2020-03-20?d=tree,map&p=grid)

* [Información básica del COVID-19](https://nextstrain.org/narratives/ncov/sit-rep/es/2020-03-20?n=2).     
* [Nota acerca del muestreo de casos](https://nextstrain.org/narratives/ncov/sit-rep/es/2020-03-20?n=3).
* [Relación entre sequenciamiento del virus y el historial de viajes](https://nextstrain.org/narratives/ncov/sit-rep/es/2020-03-20?n=6).
* [Introdución del COVID-19 en la mayoría de países del mundo](https://nextstrain.org/narratives/ncov/sit-rep/es/2020-03-20?n=8).
* [Como los brotes de una epimedia crecen y se propagan](https://nextstrain.org/narratives/ncov/sit-rep/es/2020-03-20?n=9).
* [Crecimiento del brote en el estado de Washington en EEUU](https://nextstrain.org/narratives/ncov/sit-rep/es/2020-03-20?n=10).
* [Propagación del brote en el estado de Washington en EEUU](https://nextstrain.org/narratives/ncov/sit-rep/es/2020-03-20?n=11).
* [Impacto del distaciamiento social en el número de casos](https://nextstrain.org/narratives/ncov/sit-rep/es/2020-03-20?n=12).
* [¡Qué puedes hacer tú](https://nextstrain.org/narratives/ncov/sit-rep/es/2020-03-20?n=13)!
* [Creditos Científicos](https://nextstrain.org/narratives/ncov/sit-rep/es/2020-03-20?n=14).

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown
# Resumen Ejecutivo
Este informe analiza 723 genomas públicos del COVID-19. La comparación de estos genomas virales nos permite caracterizar como el COVID-19 evoluciona y va mutando en las diferentes regiones del mundo.

En este reporte, subrayamos que el virus ha sido introducido y circula en muchas partes del mundo.
Para saber el grado de circulación local del COVID-19 así como las medidas de mitigación que pueden reducir su velocidad de transmisión se necesita conocer la dinámicas locales del brote epidémico. Esto ultimo solo puede lograrse a través de la aplicación en masa de pruebas de diagnostico.

En resumen, es fundamental tener pruebas de diagnostico rápido enfocadas en disminuir la transmisión local.

En la actualización de esta semana, reportamos que:

* Evidencia de introducciones del virus relacionadas con viajes en varias partes del mundo.
* Lugares con introducciones recientes del virus verán un incremento importante en el numero de casos en las siguientes 4-8 semanas. Deben empezar a prepararse hoy.
* Lugares que han implementado recientemente medidas de distanciamiento social seguirán viendo incrementos en el numero de casos acumulados en el corto plazo pero estos alcanzaran un máximo para luego empezar a disminuir en el mediano y largo plazo.
* Pruebas de diagnostico tanto para casos activos y en remisión/recuperados serán vitales para lidiar con la epidemia.

```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text 2 -->
# [Información básica del COVID-19](https://nextstrain.org/ncov/2020-03-20)
A continuación, hemos preparado algunos recursos de información que vale la pena leer para familiarizarse con el COVID-19 y el virus que lo causa, SARS-CoV-2. Esta información facilitará la interpretación de los datos que presentamos en este reporte.

Si no estas familiarizado con la filogenética, te recomendamos leer la siguiente introduccion ['Cómo leer philogenias'](https://nextstrain.org/narratives/trees-background/es) antes de leer el presente reporte.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

## Recursos de información relacionados con el COVID-19

<div>
  <a href="https://nextstrain.org/help/coronavirus/human-CoV"><img alt="microscopy image of coronaviruses" width="100" src="https://nextstrain.org/static/ncov_narrative-76cfd610d11ef708d213a3170de9519f.png"/> Antecendentes sobre los coronavirus </a>

  <a href="https://nextstrain.org/help/coronavirus/SARS-CoV-2"><img alt="illustration of a coronavirus" width="100" src="http://data.nextstrain.org/img_nCoV-CDC.jpg"/> Antecendentes del brote de COVID-19 </a>

  <a href="https://nextstrain.org/narratives/trees-background/es/"><img alt="cartoon of a phylogenetic tree" width="100" src="http://data.nextstrain.org/img_toy_alignment_mini.png"/> Cómo leer Filogenias </a>
</div>

## Material adicional
* Resumen del brote de SARS-CoV-2 en [Wikipedia](https://es.wikipedia.org/wiki/Pandemia_de_enfermedad_por_coronavirus_de_2019-2020).
* Los casos del COVID-19 usados en este documento son los reportados por la [OMS](https://www.who.int/docs/default-source/coronaviruse/situation-reports/20200318-sitrep-58-covid-19.pdf?sfvrsn=20876712_2) al 20-03-2020.
```

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text 3 -->
# [Nota acerca del muestreo de casos](https://nextstrain.org/ncov/2020-03-20?c=country&r=country&d=map&p=grid&legend=closed)
El presente reporte documenta la información de muestras tomadas en 36 países de 6 continentes. Este es un logro muy importante, el secuenciamiento de un virus desconocido de ARN en el medio de una pandemia, que ha sido posible gracias al trabajo sacrificado y la cooperación para compartir datos de muchos científicos y médicos en todo el mundo.
<br><br>
A pesar que estos datos nos permiten hacer importantes inferencias acerca del brote del virus y monitorear su propagación en tiempo real, debemos enfatizar que estas conclusiones son limitadas en su representatividad a nivel mundial.

Del número total de casos de COVID-19, solo una parte son diagnosticados. De los casos diagnosticados, solo una parte tienen su genoma secuenciado.
Las muestras de casos que son diagnosticados y secuenciados varían mucho entre diferentes regiones geográficas y también en el tiempo. Adicionalmente, los arboles filogenéticos de brote epidémicos tienen incertidumbre estadística intrínseca.
<br><br>
Veamos algunos ejemplos.

<!-- There is NO right-side text -->

<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text 4 -->
# [Algunas regiones estan sub-representadas en los datos actuales](https://nextstrain.org/ncov/2020-03-20?c=country&d=map&f_region=Central%20America,Oceania,South%20America,Africa&legend=closed&p=full&r=country)
El mapa muestra muy pocas secuencias del hemisferio Sur. Por ejemplo, solo se tienen secuencias de 4 de un total de 25 países de Latinoamérica donde [la OMS ha reportado](https://www.who.int/docs/default-source/coronaviruse/situation-reports/20200319-sitrep-59-covid-19.pdf?sfvrsn=c3dcdef9_2) casos de COVID-19. Esto NO se debe a que el COVID-19 no se haya propagado en otros países, o que estos casos no proporcionen información crucial para caracterizar el virus, sino simplemente que las secuencias genómicas no están disponibles en estos lugares.
<br><br>
<!-- There is NO right-side text -->


<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text 5 -->
# [Otras regiones estan sobrerepresentadas en los datos actuales](https://nextstrain.org/ncov/2020-03-20?c=country&d=map&f_region=Europe&legend=closed&p=full&r=country)
En otras áreas, como los Países Bajos, hay una gran cantidad relativa de secuencias genómicas disponibles con respecto al numero total de casos.
<br><br>
Por lo tanto, al leer este documento se debe tener en cuenta que el tamaño de los círculos en el mapa es solo proporcional a la cantidad de información disponible y no representa la magnitud del brote en cada región. Mas información sobre el efecto de esto en las inferencias estadísticas pueden encontrarse [aquí](https://nextstrain.org/narratives/trees-background/es?n=8).
<br><br>

<!-- There is NO right-side text -->

<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text 6 -->
# [Los datos de secuenciamiento genético pueden ser usados para validar los historiales de viaje](https://nextstrain.org/ncov/2020-03-20?legend=open&c=division_exposure&label=clade:A1a&d=tree)

Identificar el momento de infección de un caso es importante para entender cuales son las areas que experimentan transmisión local en comparación con aquellas donde la introducción se dio por viajes (caso cero).
<br><br>
Aquí, el color del árbol representa la historia de viaje del caso (si es conocida). Por ejemplo, la secuencia canadiense (Canada/BC_78548/2020) en el medio del árbol representa un caso con historial de viaje en Europa. Esta infección esta agrupada con otras secuencias europeas pues es casi seguro, en términos estadísticos, que es una infección relacionada con un viaje.


<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text 7 -->
# [Los datos de secuenciamiento genético pueden ser usados para validar los historiales de viaje](https://nextstrain.org/ncov/2020-03-20?c=division_exposure&d=tree&f_division_exposure=Iran&legend=open&p=full)

Si nos alejamos un poco, podemos ver que esto no siempre es tan consistente.
<br><br>
Aquí, podemos ver que casi todos los casos que han reportado algún viaje reciente a Irán forman un conglomerado en el medio del árbol.
<br><br>
Hacia la cabecera del árbol, también podemos ver que existe un caso canadiense con historial de viaje a Irán. Sin embargo, en este caso esta secuencia no se agrupa cerca del resto de casos con historial de viaje a Irán.
<br><br>
Aunque es probable que esta persona se haya infectado a través de una segunda cadena de transmisión en Irán, no podemos confirmar que esta infección se haya dado a través de un viaje hasta que no haya mas datos precisos.

<!-- There is no right side text -->

<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text 8 -->
# [Introducción del COVID-19 en la mayoría de países del mundo](https://nextstrain.org/ncov/2020-03-20?legend=closed&c=country&d=tree,map&p=grid)
En el árbol, podemos ejemplos donde los casos están inter-relacionados los unos con los otros.
Esto indica que debido al inevitable movimiento migratorio humano el virus del COVID-19 ha sido introducido en muchas partes del mundo.
<br><br>
En efecto, [la OMS reporta](https://www.who.int/docs/default-source/coronaviruse/situation-reports/20200319-sitrep-59-covid-19.pdf?sfvrsn=c3dcdef9_2) que 159 países de un total del 195 en todo el mundo tienen casos confirmados.
<br><br>
Sin embargo, no todas las introducciones resultan en epidemias.
<!-- There is no right side text -->

<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text 9 -->
# [La epidemia se hace evidente luego de 4-8 semanas desde las primeras introducciones](https://nextstrain.org/ncov/2020-03-20?legend=closed&c=num_date&d=tree&f_division=Washington&label=clade:B1&p=full)

Pareciera que el brote de COVID-19 ha explotado súbitamente esta semana.
Para mucha gente, de la noche a la mañana, la epidemia ha pasado de ser una noticia abstracta en algunos países del mundo a ser algo que impacta directamente su vida diaria.
<br><br>
Debemos recordar que el brote ha estado creciendo sin ser detectado por varias semanas. Y no se trata necesariamente que el virus haya sido introducido a lugares nuevos con mas frecuencia.
Es mas mucho probable que los países se hayan hecho conscientes de brotes locales a consecuencia de introducciones en semanas anteriores.
<br><br>
Áreas donde aun no se han detectado brotes locales deben prepararse ahora y empezar diagnósticos preliminares.

<!-- This is the right-side text -->

```auspiceMainDisplayMarkdown
# ¿Cómo las introducciones se convierten en brotes?

A veces, las introducciones del virus no resultan en casos secundarios o en brotes locales, especialmente si el caso índice (el primer caso) es rápidamente detectado y aislado. Sin embargo, el virus muchas veces se puede propagar en la comunidad local sin detección hasta que el brote llegua a un tamaño notable.

Abajo pueden ver una situación teórica que explica cómo sucedió el brote en Wuhan. El eje “y” representa el espacio geográfico; el eje “x” representa el tiempo transcurrido. El área sombreada representa el número de casos.

“COVID-19 en Wuhan se propagó desde un caso índice en noviembre del 2019 a miles de casos en mediados de enero, propagándose desde una introducción principal a un brote con transmisión local masiva en solo 10 semanas. Dado que creemos que las introducciones internacionales del virus empezaron en mediados de enero, tenemos aproximadamente 10 semanas desde entonces hasta finales de marzo para controlar los brotes individuales antes que se vuelvan masivos.” [- Trevor Bedford](https://twitter.com/trvrb/status/1226241284207038464), Feb 2020

Los brotes nuevos, a su vez, pueden ser fuentes del virus para introducciones adicionales en otros lugares del mundo.

<img src="https://github.com/nextstrain/ncov/raw/master/figures/local-spark-expansion.jpeg" width="70%">


```

<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text 10 -->
# [Cómo crecen los brotes: un ejemplo](https://nextstrain.org/ncov/2020-03-20?legend=closed&d=tree,map&f_division=Washington&label=clade:B1&p=grid&r=location)

Se puede apreciar un ejemplo claro de la propagación del virus en los datos que provienen desde el estado de Washington. Enfocándonos en la raíz del grupo principal en el árbol filogénetico, podemos concluir que el virus probablemente fue introducido a esta área entre finales de enero y mediados de febrero del 2020 ([methods](https://nextstrain.org/narratives/trees-background/es?n=6)).
<br><br>
Ahora, a mediados de marzo (~6 semanas después), podemos apreciar que el brote local en el estado de Washington está creciendo rápidamente. Usando solamente los datos de la secuenciación genómica, podemos apreciar que la tasa de duplicación del virus (lo rápido que los casos del virus se duplican) se encuentra entre 3 a 6 días, asumiendo una población con crecimiento exponencial.  
<br>
<img src="https://github.com/nextstrain/ncov/raw/master/figures/wa_doubling-rate_2020-03-19.png" width="70%">


<!-- There is no right side text -->


<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text 11-->
# [Cómo crecen los brotes: un ejemplo](https://nextstrain.org/ncov/2020-03-20?legend=closed&d=tree,map&f_division=Washington,Utah&label=clade:B1&p=grid)

También podemos apreciar un ejemplo adicional de cómo brotes que originalmente eran locales pueden servir como fuentes del virus para introducirlo en otros lugares. Las muestras aisladas de casos en el estado de Utah se encuentran en la parte superior del árbol filogénetico (en color naranja) encajadas en el grupo principal del brote de Washington. Con estos datos podemos inferir que hubo una introducción del virus desde Washington hacia Utah, tomando en cuenta que pasos intermediarios en la cadena de transmisión también son posibles.
<br><br>
No sabemos si esta introducción resultará en un brote local en Utah, pero si el virus se continúa propagando después de estas introducciones, esperamos apreciar el desarrollo del brote a través de las próximas 4 semanas.
<br><br>
Lo que presentamos es solamente un ejemplo. Es posible que hayan ocurrido otras introducciones del virus a Utah (o a otros lugares) que no aparecen en el árbol filogénetico si el virus de esos casos principales no fue secuenciado. Es altamente importante que las áreas donde los brotes locales no son todavía aparentes empiecen preparaciones y pruebas diagnósticas de vigilancia de salud pública.   

<!-- There is no right side text -->

<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text 12 -->
# [Medidas de mitigación toman tiempo, pero salvan vidas](https://nextstrain.org/ncov/2020-03-20)

Siguiendo la misma lógica anteriormente discutida, deducimos que probablemente existen muchas más cadenas de transmisión local del virus que aún no conocemos.
<br><br>
Esta conclusión significa que hasta en áreas que recientemente empezaron medidas de distanciación social, todavía esperamos que el número de casos aumente a través de las próximas semanas.

El aumento en el número de casos NO indica que estas medidas no están funcionando. Toma tiempo para que las personas que ya están infectadas (y posiblemente también los miembros de su hogar) muestren síntomas, reciban tratamiento, y se recuperen. También esperamos ver un aumento en el número de casos a medida que más pruebas diagnósticas se vuelvan más disponibles
<br><br>
Es vital mantener la distancia social durante ese tiempo. Como pueden apreciar a la derecha, el número de casos sigue aumentando justamente después de la implementación de la intervención, pero, a través del tiempo, el número total de casos disminuye considerablemente.

<!-- This is the right-side text -->

```auspiceMainDisplayMarkdown
## La distanciación social no previene todos los casos nuevos inmediatamente, pero a través del tiempo disminuye altamente el número total de casos y las muertes.
Distanciación social – es decir, disminuir el número de personas con cual te encuentras cada día – puede ser difícil, pero es increíblemente beneficioso para el público.  
 Si cada persona disminuye su número de contactos diarios por un 25%, esperamos apreciar una reducción del 50% en el número de casos acumulados en el próximo mes ([Klein et al., 2020-03-13](https://institutefordiseasemodeling.github.io/COVID-public/reports/Working%20paper%20%E2%80%93%20model-based%20estimates%20of%20COVID-19%20burden%20in%20King%20and%20Snohomish%20counties%20through%20April%207.pdf)).
<div>
  <img src="https://github.com/nextstrain/ncov/raw/master/figures/social-distancing-efficacy.png" width="70%">
</div>

```

<!-- ############ SLIDE BREAK ############# -->
<!-- This is left-side text 13-->
# [Hallazgos Principales](https://nextstrain.org/ncov/2020-03-20?c=country&d=map&p=full)
- El virus ha sido introducido a varias partes del mundo múltiples veces.  
<br>
- Podemos apreciar evidencia de transmisión local en muchos logares; más, esperamos que muchas introducciones previas del virus se desarrollen en brotes locales en las próximas semanas.  
<br>
- Controlando brotes locales a través de la distanciación social es vital para:
  - Reducir el sobrecargo de los sistemas de saludo (la campaña #aplanarlacurva #FlattenTheCurve)  
  - Disminuir el número total de casos y las muertes
  - Abonar tiempo para el desarrollo de tratamientos terapéuticos y vacunas.  

<!-- This is the right-side text -->

```auspiceMainDisplayMarkdown
# Pasos que pueden tomar los
## ...individuos
* Reducir el número de personas con cual tienes contacto diariamente, especialmente si eres parte de un grupo de gente vulnerable al virus (por ejemplo: gente mayor y gente con condiciones médicas pre-existentes).
* Recuerde que, aunque usted no sea particularmente vulnerable al virus, hay muchas personas que sí lo son; debe seguir estas prácticas para proteger a otros.
* Lávese las manos “como si acabas de picar un jalapeño y tienes que cambiar tus lentes de contacto”.  
* Quédese en casa lo más posible, especialmente si está enfermo; prepárese con provisiones adicionales por si tiene que estar en auto-cuarentena.  


## ...oficiales gubernamentales    
* Ofrecer pruebas diagnósticas gratuitas y fácilmente disponibles.  
* Establecer políticas de distanciación social.  
* Financiar e implementar seguimiento de contactos
* Apoyar financieramente a las personas y los establecimientos impactados por las políticas de distanciación social (entre ellos: trabajadores que les pagan por hora, personas responsables por gente mayor o niños pequeños, empresas pequeñas, etc.).

```

<!-- ############ SLIDE BREAK ############# -->



<!-- This is left-side text 14-->
# [Creditos Científicos](https://nextstrain.org/ncov/2020-03-20?d=map&c=author)

Nos gustaría reconocer el increíble y oportuno trabajo realizado por todos los científicos involucrados en este brote, pero particularmente aquellos que trabajan en China.
Solo mediante el intercambio rápido de datos genómicos y metadatos se pueden realizar análisis como estos.

<br>

También agradecemos a [GISAID](https://gisaid.org) por proporcionar la plataforma a través de la cual estos datos se pueden depositar y compartir.

<!-- Do not need to translate institutions names -->
<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

Agradecemos los datos generados por estos laboratorios:

* Arizona Department of Health Services
* Auckland Hospital
* BCCDC Public Health Laboratory
* Bamrasnaradura Hospital
* Beijing Institute of Microbiology and Epidemiology
* Bundeswehr Institute of Microbiology
* CNR Virus des Infections Respiratoires - France SUD
* CR&WISCO GENERAL HOSPITAL
* California Department of Health
* California Department of Public Health
* Center of Medical Microbiology, Virology, and Hospital Hygiene
* Center of Medical Microbiology, Virology, and Hospital Hygiene, University of Duesseldorf
* Centers for Disease Control, R.O.C. (Taiwan)
* Centre for Human and Zoonotic Virology (CHAZVY), College of Medicine University of Lagos/Lagos University Teaching Hospital (LUTH), part of the Laboratory Network of the Nigeria Centre for Disease Control (NCDC)
* Centre for Infectious Diseases and Microbiology - Public Health
* Centre for Infectious Diseases and Microbiology Laboratory Services
* Centre for Infectious Diseases and Microbiology- Public Health
* Centro Hospital do Porto, E.P.E. - H. Geral de Santo Antonio
* Centro Hospitalar e Universitario de Sao Joao, Porto
* Charite Universitatsmedizin Berlin, Institute of Virology; Institut fur Mikrobiologie der Bundeswehr, Munich
* Department of Infectious Diseases, Istituto Superiore di Sanita, Roma , Italy
* Department of Infectious Diseases, Istituto Superiore di Sanita, Rome, Italy
* Department of Infectious and Tropical Diseases, Bichat Claude Bernard Hospital, Paris
* Department of Internal Medicine, Triemli Hospital
* Department of Laboratory Medicine, National Taiwan University Hospital
* Department of Microbiology, Institute for Viral Diseases, College of Medicine, Korea University
* Department of Pathology, Toshima Hospital
* Department of Virology III, National Institute of Infectious Diseases
* Department of Virology and Immunology, University of Helsinki and Helsinki University Hospital, Huslab Finland
* Department of microbiology laboratory,Anhui Provincial Center for Disease Control and Prevention
* Dept. of Pathology, National Institute of Infectious Diseases
* Dept. of Virology III, National Institute of Infectious Diseases
* Dienst Gezondheid & Jeugd Zuid-Holland Zuid
* Division of Infectious Diseases, Department of Internal Medicine, Korea University College of Medicine
* Division of Infectious Diseases, University Hospital Zurich
* Division of Viral Diseases, Center for Laboratory Control of Infectious Diseases, Korea Centers for Diseases Control and Prevention
* Dutch COVID-19 response team
* ErasmusMC
* Foundation Elisabeth-Tweesteden Ziekenhuis
* Foundation Pamm
* Fujian Center for Disease Control and Prevention
* General Hospital of Central Theater Command of People's Liberation Army of China
* Guangdong Provincial Center for Diseases Control and Prevention; Guangdong Provincial Public Health
* Guangdong Provincial Center for Diseases Control and Prevention; Guangdong Provinical Public Health
* Guangdong Provincial Center for Diseases Control and Prevention;Guangdong Provincial Institute of Public Health
* Guangdong Provincial Institution of Public Health, Guangdong Provinical Center for Disease Control and Prevention
* HUS Diagnostiikkakeskus, Hallinto
* Hangzhou Center for Disease Control and Prevention
* Hangzhou Center for Disease and Control Microbiology Lab
* Harborview Medical Center
* Hong Kong Department of Health
* Hospital Israelita Albert Einstein
* Hospital Sao Joaquim Beneficencia Portuguesa
* IL Department of Public Health Chicago Laboratory
* INMI Lazzaro Spallanzani IRCCS
* Indian Council of Medical Research - National Institute of Virology
* Indian Council of Medical Research-National Institute of Virology
* Institute of Pathogen Biology, Chinese Academy of Medical Sciences & Peking Union Medical College
* Institute of Viral Disease Control and Prevention, China CDC
* Instituto Nacional de Enfermedades Respiratorias
* Jingzhou Center for Disease Control and Prevention
* KU Leuven, Clinical and Epidemiological Virology
* Klinik Hirslanden Zurich
* Korea Centers for Disease Control & Prevention (KCDC) Center for Laboratory Control of Infectious Diseases Division of Viral Diseases
* Laboratoire National de Sante
* Laboratoire de Virologie, HUG
* Laboratorio di Microbiologia e Virologia, Universita Vita-Salute San Raffaele, Milano
* Laboratory Medicine
* Lapland Central Hospital
* MHC Brabant Zuidoost
* MHC Drente
* MHC Flevoland
* MHC Gooi & Vechtstreek
* MHC Haaglanden
* MHC Hart voor Brabant
* MHC Kennemerland
* MHC Rotterdam-Rijnmond
* MHC Utrecht
* MHC West-Brabant
* MSHS Clinical Microbiology Laboratories
* Massachusetts Department of Public Health
* Monash Medical Centre
* NHC Key laboratory of Enteric Pathogenic Microbiology, Institute of Pathogenic Microbiology
* National Centre for Infectious Diseases
* National Influenza Center - National Institute of Hygiene and Epidemiology (NIHE)
* National Influenza Centre, National Public Health Laboratory, Kathmandu, Nepal
* National Institute for Viral Disease Control and Prevention, China CDC
* National Public Health Laboratory
* National Public Health Laboratory, National Centre for Infectious Diseases
* Pathology Queensland
* Providence Regional Medical Center
* Public Health Ontario Laboratory
* RIVM
* Respiratory Virus Unit, Microbiology Services Colindale, Public Health England
* Seattle Flu Study
* Second Hospital of Anhui Medical University
* Serology, Virology and OTDS Laboratories (SAViD), NSW Health Pathology Randwick
* Servicio Microbiologia. Hospital Clinico Universitario. Valencia.
* Shenzhen Key Laboratory of Pathogen and Immunity, National Clinical Research Center for Infectious Disease, Shenzhen Third People's Hospital
* Shenzhen Third People's Hospital
* Singapore General Hospital
* Singapore General Hospital, Molecular Laboratory, Division of Pathology
* Sorbonne Universite, Inserm et Assistance Publique-Hopitaux de Paris (Pitie Salpetriere)
* South China Agricultural University
* State Health Office Baden-Wuerttemberg
* Taiwan Centers for Disease Control
* Texas Department of State Health Services
* The Central Hospital Of Wuhan
* The National Institute of Public Health Center for Epidemiology and Microbiology
* The University of Hong Kong - Shenzhen Hospital
* Tianmen Center for Disease Control and Prevention
* UCD National Virus Reference Laboratory
* UW Virology Lab
* Union Hospital of Tongji Medical College, Huazhong University of Science and Technology
* Unknown
* Valley Medical Center
* Virology Department, Sheffield Teaching Hospitals NHS Foundation Trust
* Virology Unit, Institut Pasteur du Cambodge.
* WA State Department of Health
* Wales Specialist Virology Centre
* Washington State Department of Health
* Washington State Public Health Lab
* Weifang Center for Disease Control and Prevention
* West of Scotland Specialist Virology Centre, NHSGGC
* Wisconsin Department of Health Services
* Wuhan Fourth Hospital
* Wuhan Institute of Virology, Chinese Academy of Sciences
* Wuhan Jinyintan Hospital
* Wuhan Lung Hospital
* Yongchuan District Center for Disease Control and Prevention
* Zhejiang Provincial Center for Disease Control and Prevention
* Zhongxian Center for Disease Control and Prevention
* Andersen Lab, The Scripps Research Institute
* Arizona Department of Health Services
* Auckland Hospital
* BCCDC Public Health Laboratory
* Bamrasnaradura Hospital
* Beijing Institute of Microbiology and Epidemiology
* Bundeswehr Institute of Microbiology
* CNR Virus des Infections Respiratoires - France SUD
* CR&WISCO GENERAL HOSPITAL
* California Department of Health
* California Department of Public Health
* Center of Medical Microbiology, Virology, and Hospital Hygiene, University of Duesseldorf
* Centers for Disease Control, R.O.C. (Taiwan)
* Centre Hositalier Universitaire de Rouen Laboratoire de Virologie
* Centre Hospitalier Compiegne Laboratoire de Biologie
* Centre Hospitalier Regional Universitaire de Nantes Laboratoire de Virologie
* Centre Hospitalier Rene Dubois Laboratoire de Microbiologie - Bat A
* Centre for Human and Zoonotic Virology (CHAZVY), College of Medicine University of Lagos/Lagos University Teaching Hospital (LUTH), part of the Laboratory Network of the Nigeria Centre for Disease Control (NCDC)
* Centre for Infectious Diseases and Microbiology - Public Health
* Centre for Infectious Diseases and Microbiology Laboratory Services
* Centre for Infectious Diseases and Microbiology- Public Health
* Centro Hospital do Porto, E.P.E. - H. Geral de Santo Antonio
* Centro Hospitalar e Universitario de Sao Joao, Porto
* Charite Universitatsmedizin Berlin, Institute of Virology; Institut fur Mikrobiologie der Bundeswehr, Munich
* Clinica Alemana de Santiago, Chile
* Clinica Santa Maria, Santiago, Chile
* Department of Infectious Diseases, Istituto Superiore di Sanita, Roma , Italy
* Department of Infectious Diseases, Istituto Superiore di Sanita, Rome, Italy
* Department of Infectious and Tropical Diseases, Bichat Claude Bernard Hospital, Paris
* Department of Internal Medicine, Triemli Hospital
* Department of Laboratory Medicine, National Taiwan University Hospital
* Department of Microbiology, Institute for Viral Diseases, College of Medicine, Korea University
* Department of Pathology, Toshima Hospital
* Department of Virology III, National Institute of Infectious Diseases
* Department of Virology and Immunology, University of Helsinki and Helsinki University Hospital, Huslab Finland
* Department of Virus and Microbiological Special diagnostics, Statens Serum Institut, Copenhagen, Denmark.
* Department of microbiology laboratory,Anhui Provincial Center for Disease Control and Prevention
* Dept. of Pathology, National Institute of Infectious Diseases
* Dept. of Virology III, National Institute of Infectious Diseases
* Dienst Gezondheid & Jeugd Zuid-Holland Zuid
* Division of Infectious Diseases, Department of Internal Medicine, Korea University College of Medicine
* Division of Infectious Diseases, University Hospital Zurich
* Division of Viral Diseases, Center for Laboratory Control of Infectious Diseases, Korea Centers for Diseases Control and Prevention
* Dutch COVID-19 response team
* ErasmusMC
* Foundation Elisabeth-Tweesteden Ziekenhuis
* Foundation Pamm
* Fujian Center for Disease Control and Prevention
* General Hospital of Central Theater Command of People's Liberation Army of China
* Gorgas Memorial Institute for Health Studies
* Guangdong Provincial Center for Diseases Control and Prevention; Guangdong Provincial Public Health
* Guangdong Provincial Center for Diseases Control and Prevention; Guangdong Provinical Public Health
* Guangdong Provincial Center for Diseases Control and Prevention;Guangdong Provincial Institute of Public Health
* Guangdong Provincial Institution of Public Health, Guangdong Provinical Center for Disease Control and Prevention
* HUS Diagnostiikkakeskus, Hallinto
* Hangzhou Center for Disease Control and Prevention
* Hangzhou Center for Disease and Control Microbiology Lab
* Harborview Medical Center
* Hong Kong Department of Health
* Hopital Instruction des Armees - BEGIN
* Hopital Robert Debre Laboratoire de Virologie
* Hopitaux universitaires de Geneve Laboratoire de Virologie
* Hospital Israelita Albert Einstein
* Hospital Sao Joaquim Beneficencia Portuguesa
* Hospital de Talca, Chile
* IL Department of Public Health Chicago Laboratory
* INMI Lazzaro Spallanzani IRCCS
* Indian Council of Medical Research - National Institute of Virology
* Indian Council of Medical Research-National Institute of Virology
* Institute of Pathogen Biology, Chinese Academy of Medical Sciences & Peking Union Medical College
* Institute of Viral Disease Control and Prevention, China CDC
* Instituto Nacional de Enfermedades Respiratorias
* Jingzhou Center for Disease Control and Prevention
* KU Leuven, Clincal and Epidemiological Virology
* KU Leuven, Clinical and Epidemiological Virology
* Klinik Hirslanden Zurich
* Korea Centers for Disease Control & Prevention (KCDC) Center for Laboratory Control of Infectious Diseases Division of Viral Diseases
* LACEN RJ - Laboratorio Central de Saude Publica Noel Nutels
* LACEN/ES - Laboratorio Central de Saude Publica do Espirito Santo
* Laboratoire National de Sante
* Laboratoire de Virologie Institut de Virologie - INSERM U 1109 Hopitaux Universitaires de Strasbourg
* Laboratoire de Virologie, HUG
* Laboratorio Central de Saude Publica Professor Goncalo Moniz  LACEN/BA
* Laboratorio di Microbiologia e Virologia, Universita Vita-Salute San Raffaele, Milano
* Laboratory Medicine
* Laboratory of Molecular Virology, Pontificia Universidad Catolica de Chile
* Lapland Central Hospital
* MHC Brabant Zuidoost
* MHC Drente
* MHC Flevoland
* MHC Gooi & Vechtstreek
* MHC Haaglanden
* MHC Hart voor Brabant
* MHC Kennemerland
* MHC Rotterdam-Rijnmond
* MHC Utrecht
* MHC West-Brabant
* MSHS Clinical Microbiology Laboratories
* Massachusetts Department of Public Health
* Minnesota Department of Health, Public Health Laboratory
* Monash Medical Centre
* NHC Key laboratory of Enteric Pathogenic Microbiology, Institute of Pathogenic Microbiology
* NYU Langone Health
* National Centre for Infectious Diseases
* National Influenza Center - National Institute of Hygiene and Epidemiology (NIHE)
* National Influenza Centre, National Public Health Laboratory, Kathmandu, Nepal
* National Institute for Viral Disease Control and Prevention, China CDC
* National Public Health Laboratory
* National Public Health Laboratory, National Centre for Infectious Diseases
* Pathology Queensland
* Providence Regional Medical Center
* Public Health Ontario Laboratory
* R. G. Lugar Center for Public Health Research,  National Center for Disease Control and Public Health (NCDC) of Georgia.
* RIVM
* Regional Virus Laboratory, Belfast
* Respiratory Virus Unit, Microbiology Services Colindale, Public Health England
* Seattle Flu Study
* Second Hospital of Anhui Medical University
* Serology, Virology and OTDS Laboratories (SAViD), NSW Health Pathology Randwick
* Servicio Microbiologia, Hospital Clinico Universitario, Valencia
* Servicio Microbiologia. Hospital Clinico Universitario. Valencia.
* Shandong Provincial Center for Disease Control and Prevention
* Shenzhen Key Laboratory of Pathogen and Immunity, National Clinical Research Center for Infectious Disease, Shenzhen Third People's Hospital
* Shenzhen Third People's Hospital
* Singapore General Hospital
* Singapore General Hospital, Molecular Laboratory, Division of Pathology
* Sorbonne Universite, Inserm et Assistance Publique-Hopitaux de Paris (Pitie Salpetriere)
* South China Agricultural University
* State Health Office Baden-Wuerttemberg
* State Key Laboratory for Diagnosis and Treatment of Infectious Diseases, National Clinical Research Center for Infectious Diseases, First Affiliated Hospital, Zhejiang University School of Medicine, Hangzhou, China. 310003
* State Key Laboratory of Respiratory Disease, National Clinical Research Center for Respiratory Disease, Guangzhou Institute of Respiratory Health, the First Affiliated Hospital of Guangzhou Medical University
* Tai Lung Veterinary Laboratory, Agriculture, Fisheries and Conservation Department
* Taiwan Centers for Disease Control
* Texas Department of State Health Services
* The Central Hospital Of Wuhan
* The National Institute of Public Health Center for Epidemiology and Microbiology
* The University of Hong Kong - Shenzhen Hospital
* Tianmen Center for Disease Control and Prevention
* UCD National Virus Reference Laboratory
* UW Virology Lab
* Union Hospital of Tongji Medical College, Huazhong University of Science and Technology
* Unknown
* Utah Public Health Laboratory
* Valley Medical Center
* Viral Respiratory Lab, National Institute for Biomedical Research (INRB)
* Virology Department, Royal Infirmary of Edinburgh, NHS Lothian
* Virology Department, Sheffield Teaching Hospitals NHS Foundation Trust
* Virology Unit, Institut Pasteur du Cambodge.
* WA State Department of Health
* WHO National Influenza Centre Russian Federation
* Wales Specialist Virology Centre
* Washington State Department of Health
* Washington State Public Health Lab
* Weifang Center for Disease Control and Prevention
* West of Scotland Specialist Virology Centre, NHSGGC
* Wisconsin Department of Health Services
* Wuhan Fourth Hospital
* Wuhan Institute of Virology, Chinese Academy of Sciences
* Wuhan Jinyintan Hospital
* Wuhan Lung Hospital
* Yongchuan District Center for Disease Control and Prevention
* Zhejiang Provincial Center for Disease Control and Prevention
* Zhongxian Center for Disease Control and Prevention

```


<!-- ############ SLIDE BREAK ############# -->

<!-- This is left-side text 13-->
# [Créditos científicos detallados](https://nextstrain.org/ncov/2020-03-20?d=map&c=author)

Los datos utilizados en este reporte son datos públicos compartidos por cada laboratorio a través de [GISAID](https://gisaid.org).

Agradecemos cordialmente sus contribuciones.
<br>

En el panel de la derecha especifícamos las secuencias compartidas por cada laboratorio.

<!-- This is right-side text -->
```auspiceMainDisplayMarkdown

Los genomas del SARS-CoV-2 fueron generosamente compartidos por los científicos en estos laboratorios:

* Arizona Department of Health Services
	* USA/AZ1/2020

* Auckland Hospital
	* NewZealand/01/2020

* BCCDC Public Health Laboratory
	* Canada/BC_37_0-2/2020

* Bamrasnaradura Hospital
	* Nonthaburi/61/2020
	* Nonthaburi/74/2020

* Beijing Institute of Microbiology and Epidemiology
	* pangolin/Guangdong/P2S/2019
	* pangolin/Guangxi/P1E/2017
	* pangolin/Guangxi/P2V/2017
	* pangolin/Guangxi/P3B/2017
	* pangolin/Guangxi/P4L/2017
	* pangolin/Guangxi/P5E/2017
	* pangolin/Guangxi/P5L/2017

* Bundeswehr Institute of Microbiology
	* Germany/BavPat2/2020
	* Germany/BavPat3/2020

* CNR Virus des Infections Respiratoires - France SUD
	* France/RA739/2020

* CR&WISCO GENERAL HOSPITAL
	* Wuhan/HBCDC-HB-05/2020

* California Department of Health
	* USA/CA3/2020
	* USA/CA4/2020
	* USA/CA5/2020

* California Department of Public Health
	* USA/CA-CDPH-UC1/2020
	* USA/CA-CDPH-UC2/2020
	* USA/CA-CDPH-UC3/2020
	* USA/CA-CDPH-UC4/2020
	* USA/CA-CDPH-UC5/2020
	* USA/CA-CDPH-UC6/2020
	* USA/CA-CDPH-UC7/2020
	* USA/CA-CDPH-UC8/2020
	* USA/CA-CDPH-UC9/2020
	* USA/CA1/2020
	* USA/CA2/2020
	* USA/CA6/2020
	* USA/CA7/2020
	* USA/CA8/2020
	* USA/CA9/2020
	* USA/UC-CDPH-UC11/2020

* Center of Medical Microbiology, Virology, and Hospital Hygiene
	* Germany/NRW-01/2020
	* Germany/NRW-02-1/2020
	* Germany/NRW-03/2020
	* Germany/NRW-04/2020

* Center of Medical Microbiology, Virology, and Hospital Hygiene, University of Duesseldorf
	* Germany/NRW-011/2020
	* Germany/NRW-05/2020
	* Germany/NRW-06/2020
	* Germany/NRW-07/2020
	* Germany/NRW-08/2020
	* Germany/NRW-09/2020
	* Germany/NRW-10/2020

* Centers for Disease Control, R.O.C. (Taiwan)
	* Taiwan/2/2020

* Centre for Human and Zoonotic Virology (CHAZVY), College of Medicine University of Lagos/Lagos University Teaching Hospital (LUTH), part of the Laboratory Network of the Nigeria Centre for Disease Control (NCDC)
	* Nigeria/Lagos01/2020

* Centre for Infectious Diseases and Microbiology - Public Health
	* Australia/NSW10/2020
	* Australia/NSW12/2020
	* Australia/NSW13/2020
	* Australia/NSW14/2020

* Centre for Infectious Diseases and Microbiology Laboratory Services
	* Australia/NSW01/2020
	* Australia/NSW05/2020
	* Australia/NSW06/2020
	* Australia/NSW07/2020
	* Australia/NSW08/2020
	* Australia/NSW09/2020
	* Sydney/2/2020

* Centre for Infectious Diseases and Microbiology- Public Health
	* Australia/NSW11/2020

* Centro Hospital do Porto, E.P.E. - H. Geral de Santo Antonio
	* Portugal/CV62/2020

* Centro Hospitalar e Universitario de Sao Joao, Porto
	* Portugal/CV63/2020

* Charite Universitatsmedizin Berlin, Institute of Virology; Institut fur Mikrobiologie der Bundeswehr, Munich
	* Germany/BavPat1/2020

* Department of Infectious Diseases, Istituto Superiore di Sanita, Roma , Italy
	* Italy/CDG1/2020

* Department of Infectious Diseases, Istituto Superiore di Sanita, Rome, Italy
	* Italy/SPL1/2020

* Department of Infectious and Tropical Diseases, Bichat Claude Bernard Hospital, Paris
	* France/IDF0372-isl/2020
	* France/IDF0372/2020
	* France/IDF0373/2020
	* France/IDF0386-islP1/2020
	* France/IDF0386-islP3/2020
	* France/IDF0515-isl/2020
	* France/IDF0515/2020
	* France/IDF0571/2020

* Department of Internal Medicine, Triemli Hospital
	* Switzerland/1000477102/2020
	* Switzerland/1000477377/2020

* Department of Laboratory Medicine, National Taiwan University Hospital
	* Taiwan/NTU01/2020
	* Taiwan/NTU02/2020
	* Taiwan/NTU03/2020

* Department of Microbiology, Institute for Viral Diseases, College of Medicine, Korea University
	* SouthKorea/KUMC01/2020
	* SouthKorea/KUMC02/2020
	* SouthKorea/KUMC04/2020
	* SouthKorea/KUMC06/2020

* Department of Pathology, Toshima Hospital
	* Japan/TK/20-31-3/2020

* Department of Virology III, National Institute of Infectious Diseases
	* Japan/AI/I-004/2020

* Department of Virology and Immunology, University of Helsinki and Helsinki University Hospital, Huslab Finland
	* Finland/FIN01032020/2020
	* Finland/FIN03032020A/2020
	* Finland/FIN03032020B/2020
	* Finland/FIN03032020C/2020

* Department of microbiology laboratory,Anhui Provincial Center for Disease Control and Prevention
	* Anhui/SZ005/2020

* Dept. of Pathology, National Institute of Infectious Diseases
	* Japan/NA-20-05-1/2020
	* Japan/OS-20-07-1/2020

* Dept. of Virology III, National Institute of Infectious Diseases
	* Japan/KY-V-029/2020
	* Japan/TY-WK-012/2020
	* Japan/TY-WK-501/2020
	* Japan/TY-WK-521/2020

* Dienst Gezondheid & Jeugd Zuid-Holland Zuid
	* Netherlands/Hardinxveld_Giessendam_1364806/2020

* Division of Infectious Diseases, Department of Internal Medicine, Korea University College of Medicine
	* SouthKorea/KUMC03/2020
	* SouthKorea/KUMC05/2020

* Division of Infectious Diseases, University Hospital Zurich
	* Switzerland/1000477796/2020
	* Switzerland/1000477797/2020
	* Switzerland/1000477806/2020

* Division of Viral Diseases, Center for Laboratory Control of Infectious Diseases, Korea Centers for Diseases Control and Prevention
	* SouthKorea/KCDC05/2020
	* SouthKorea/KCDC06/2020
	* SouthKorea/KCDC07/2020
	* SouthKorea/KCDC12/2020
	* SouthKorea/KCDC24/2020

* Dutch COVID-19 response team
	* Netherlands/Gelderland_1/2020
	* Netherlands/Limburg_2/2020
	* Netherlands/Limburg_3/2020
	* Netherlands/Limburg_4/2020
	* Netherlands/Limburg_5/2020
	* Netherlands/Limburg_6/2020
	* Netherlands/NoordBrabant_1/2020
	* Netherlands/NoordBrabant_10/2020
	* Netherlands/NoordBrabant_11/2020
	* Netherlands/NoordBrabant_12/2020
	* Netherlands/NoordBrabant_13/2020
	* Netherlands/NoordBrabant_14/2020
	* Netherlands/NoordBrabant_15/2020
	* Netherlands/NoordBrabant_16/2020
	* Netherlands/NoordBrabant_17/2020
	* Netherlands/NoordBrabant_18/2020
	* Netherlands/NoordBrabant_19/2020
	* Netherlands/NoordBrabant_2/2020
	* Netherlands/NoordBrabant_20/2020
	* Netherlands/NoordBrabant_21/2020
	* Netherlands/NoordBrabant_22/2020
	* Netherlands/NoordBrabant_23/2020
	* Netherlands/NoordBrabant_24/2020
	* Netherlands/NoordBrabant_25/2020
	* Netherlands/NoordBrabant_26/2020
	* Netherlands/NoordBrabant_27/2020
	* Netherlands/NoordBrabant_28/2020
	* Netherlands/NoordBrabant_29/2020
	* Netherlands/NoordBrabant_3/2020
	* Netherlands/NoordBrabant_30/2020
	* Netherlands/NoordBrabant_31/2020
	* Netherlands/NoordBrabant_32/2020
	* Netherlands/NoordBrabant_33/2020
	* Netherlands/NoordBrabant_34/2020
	* Netherlands/NoordBrabant_35/2020
	* Netherlands/NoordBrabant_36/2020
	* Netherlands/NoordBrabant_37/2020
	* Netherlands/NoordBrabant_38/2020
	* Netherlands/NoordBrabant_39/2020
	* Netherlands/NoordBrabant_4/2020
	* Netherlands/NoordBrabant_5/2020
	* Netherlands/NoordBrabant_6/2020
	* Netherlands/NoordHolland_1/2020
	* Netherlands/NoordHolland_2/2020
	* Netherlands/Overijssel_1/2020
	* Netherlands/Overijssel_2/2020
	* Netherlands/Utrecht_1/2020
	* Netherlands/Utrecht_10/2020
	* Netherlands/Utrecht_11/2020
	* Netherlands/Utrecht_12/2020
	* Netherlands/Utrecht_13/2020
	* Netherlands/Utrecht_14/2020
	* Netherlands/Utrecht_15/2020
	* Netherlands/Utrecht_16/2020
	* Netherlands/Utrecht_2/2020
	* Netherlands/Utrecht_3/2020
	* Netherlands/Utrecht_4/2020
	* Netherlands/Utrecht_5/2020
	* Netherlands/Utrecht_6/2020
	* Netherlands/Utrecht_7/2020
	* Netherlands/Utrecht_8/2020
	* Netherlands/ZuidHolland_1/2020
	* Netherlands/ZuidHolland_10/2020
	* Netherlands/ZuidHolland_11/2020
	* Netherlands/ZuidHolland_13/2020
	* Netherlands/ZuidHolland_14/2020
	* Netherlands/ZuidHolland_15/2020
	* Netherlands/ZuidHolland_16/2020
	* Netherlands/ZuidHolland_17/2020
	* Netherlands/ZuidHolland_18/2020
	* Netherlands/ZuidHolland_19/2020
	* Netherlands/ZuidHolland_2/2020
	* Netherlands/ZuidHolland_20/2020
	* Netherlands/ZuidHolland_21/2020
	* Netherlands/ZuidHolland_22/2020
	* Netherlands/ZuidHolland_23/2020
	* Netherlands/ZuidHolland_24/2020
	* Netherlands/ZuidHolland_5/2020
	* Netherlands/ZuidHolland_6/2020
	* Netherlands/ZuidHolland_7/2020
	* Netherlands/ZuidHolland_8/2020
	* Netherlands/ZuidHolland_9/2020

* ErasmusMC
	* Netherlands/Nieuwendijk_1363582/2020
	* Netherlands/Rotterdam_1363790/2020

* Foundation Elisabeth-Tweesteden Ziekenhuis
	* Netherlands/Tilburg_1363354/2020
	* Netherlands/Tilburg_1364286/2020

* Foundation Pamm
	* Netherlands/Berlicum_1363564/2020

* Fujian Center for Disease Control and Prevention
	* Fujian/13/2020
	* Fujian/8/2020

* General Hospital of Central Theater Command of People's Liberation Army of China
	* Wuhan/WH01/2019
	* Wuhan/WH02/2019
	* Wuhan/WH03/2020
	* Wuhan/WH04/2020

* Guangdong Provincial Center for Diseases Control and Prevention; Guangdong Provincial Public Health
	* Foshan/20SF207/2020
	* Foshan/20SF210/2020
	* Foshan/20SF211/2020
	* Guangdong/20SF012/2020
	* Guangdong/20SF013/2020
	* Guangdong/20SF014/2020
	* Guangdong/20SF025/2020
	* Guangdong/20SF028/2020
	* Guangdong/20SF040/2020

* Guangdong Provincial Center for Diseases Control and Prevention; Guangdong Provinical Public Health
	* Guangdong/20SF174/2020
	* Guangzhou/20SF206/2020

* Guangdong Provincial Center for Diseases Control and Prevention;Guangdong Provincial Institute of Public Health
	* Guangdong/20SF201/2020

* Guangdong Provincial Institution of Public Health, Guangdong Provinical Center for Disease Control and Prevention
	* Guangdong/2020XN4239-P0034/2020
	* Guangdong/2020XN4243-P0035/2020
	* Guangdong/2020XN4273-P0036/2020
	* Guangdong/2020XN4276-P0037/2020
	* Guangdong/2020XN4291-P0038/2020
	* Guangdong/2020XN4373-P0039/2020
	* Guangdong/2020XN4433-P0040/2020
	* Guangdong/2020XN4448-P0002/2020
	* Guangdong/2020XN4459-P0041/2020
	* Guangdong/2020XN4475-P0042/2020
	* Guangdong/DG-S2-P0054/2020
	* Guangdong/DG-S41-P0056/2020
	* Guangdong/DG-S6-P0055/2020
	* Guangdong/DG-S9-P0045/2020
	* Guangdong/FS-S29-P0051/2020
	* Guangdong/FS-S30-P0052/2020
	* Guangdong/FS-S34-P0015/2020
	* Guangdong/FS-S42-P0046/2020
	* Guangdong/FS-S48-P0047/2020
	* Guangdong/FS-S50-P0053/2020
	* Guangdong/GD2020012-P0022/2020
	* Guangdong/GD2020016-P0011/2020
	* Guangdong/GD2020080-P0010/2020
	* Guangdong/GD2020085-P0043/2020
	* Guangdong/GD2020086-P0021/2020
	* Guangdong/GD2020087-P0008/2020
	* Guangdong/GD2020115-P0009/2020
	* Guangdong/GD2020134-P0031/2020
	* Guangdong/GD2020139-P0007/2020
	* Guangdong/GD2020227-P0029/2020
	* Guangdong/GD2020233-P0027/2020
	* Guangdong/GD2020234-P0023/2020
	* Guangdong/GD2020241-P0013/2020
	* Guangdong/GD2020246-P0028/2020
	* Guangdong/GD2020258-P0018/2020
	* Guangdong/GDFS2020052-P0025/2020
	* Guangdong/GDFS2020054-P0005/2020
	* Guangdong/GDFS2020056-P0044/2020
	* Guangdong/GDFS2020127-P0026/2020
	* Guangdong/GDSZ202004-P0004/2020
	* Guangdong/GDSZ202008-P0020/2020
	* Guangdong/GDSZ202009-P0032/2020
	* Guangdong/GDSZ202013-P0014/2020
	* Guangdong/GDSZ202015-P0019/2020
	* Guangdong/GZ-S6-P0050/2020
	* Guangdong/JM-S1-P0062/2020
	* Guangdong/MM-S1-P0048/2020
	* Guangdong/SZ-N128-P0057/2020
	* Guangdong/SZ-N59-P0049/2020
	* Guangdong/ZH-N22-P0059/2020
	* Guangdong/ZH-S33-P0058/2020
	* Guangdong/ZQ-S2-P0061/2020
	* Guangdong/ZS-S6-P0060/2020

* HUS Diagnostiikkakeskus, Hallinto
	* Finland/FIN-25/2020

* Hangzhou Center for Disease Control and Prevention
	* Hangzhou/HZCDC0001/2020

* Hangzhou Center for Disease and Control Microbiology Lab
	* Hangzhou/HZ-1/2020

* Harborview Medical Center
	* USA/WA3-UW1/2020
	* USA/WA9-UW6/2020

* Hong Kong Department of Health
	* HongKong/VB20024950/2020
	* HongKong/VB20026565/2020
	* HongKong/VM20001061/2020
	* HongKong/case42_VM20002493/2020
	* HongKong/case48_VM20002507/2020
	* HongKong/case52_VM20002582/2020
	* HongKong/case78_VM20002849/2020
	* HongKong/case85_VM20002868/2020
	* HongKong/case90_VM20002907/2020
	* canine/HongKong/20-02756/2020

* Hospital Israelita Albert Einstein
	* Brazil/SPBR-01/2020
	* Brazil/SPBR-02/2020
	* Brazil/SPBR-03/2020

* Hospital Sao Joaquim Beneficencia Portuguesa
	* Brazil/SPBR-04/2020
	* Brazil/SPBR-05/2020
	* Brazil/SPBR-06/2020

* IL Department of Public Health Chicago Laboratory
	* USA/IL1/2020
	* USA/IL2/2020

* INMI Lazzaro Spallanzani IRCCS
	* Italy/INMI1-cs/2020
	* Italy/INMI1-isl/2020

* Indian Council of Medical Research - National Institute of Virology
	* India/1-27/2020

* Indian Council of Medical Research-National Institute of Virology
	* India/1-31/2020

* Institute of Pathogen Biology, Chinese Academy of Medical Sciences & Peking Union Medical College
	* Wuhan/IPBCAMS-WH-01/2019
	* Wuhan/IPBCAMS-WH-02/2019
	* Wuhan/IPBCAMS-WH-03/2019
	* Wuhan/IPBCAMS-WH-04/2019
	* Wuhan/IPBCAMS-WH-05/2020

* Institute of Viral Disease Control and Prevention, China CDC
	* Wuhan/IVDC-HB-envF13-20/2020
	* Wuhan/IVDC-HB-envF13-21/2020
	* Wuhan/IVDC-HB-envF13/2020
	* Wuhan/IVDC-HB-envF54/2020

* Instituto Nacional de Enfermedades Respiratorias
	* Mexico/CDMX/InDRE_01/2020

* Jingzhou Center for Disease Control and Prevention
	* Jingzhou/HBCDC-HB-01/2020

* KU Leuven, Clinical and Epidemiological Virology
	* Belgium/GHB-03021/2020

* Klinik Hirslanden Zurich
	* Switzerland/1000477757/2020

* Korea Centers for Disease Control & Prevention (KCDC) Center for Laboratory Control of Infectious Diseases Division of Viral Diseases
	* SouthKorea/KCDC03/2020

* Laboratoire National de Sante
	* Luxembourg/Lux1/2020

* Laboratoire de Virologie, HUG
	* Switzerland/AG0361/2020
	* Switzerland/BL0902/2020
	* Switzerland/GE3121/2020
	* Switzerland/GE3895/2020
	* Switzerland/GE5373/2020
	* Switzerland/GE9586/2020
	* Switzerland/TI9486/2020
	* Switzerland/VD5615/2020

* Laboratorio di Microbiologia e Virologia, Universita Vita-Salute San Raffaele, Milano
	* Italy/UniSR1/2020

* Laboratory Medicine
	* Taiwan/CGMH-CGU-01/2020

* Lapland Central Hospital
	* Finland/1/2020

* MHC Brabant Zuidoost
	* Netherlands/Eindhoven_1363782/2020

* MHC Drente
	* Netherlands/Dalen_1363624/2020

* MHC Flevoland
	* Netherlands/Zeewolde_1365080/2020

* MHC Gooi & Vechtstreek
	* Netherlands/Blaricum_1364780/2020
	* Netherlands/Naarden_1364774/2020

* MHC Haaglanden
	* Netherlands/Nootdorp_1364222/2020

* MHC Hart voor Brabant
	* Netherlands/Oisterwijk_1364072/2020

* MHC Kennemerland
	* Netherlands/Haarlem_1363688/2020

* MHC Rotterdam-Rijnmond
	* Netherlands/Rotterdam_1364040/2020

* MHC Utrecht
	* Netherlands/Utrecht_1363564/2020
	* Netherlands/Utrecht_1363628/2020
	* Netherlands/Utrecht_1364066/2020

* MHC West-Brabant
	* Netherlands/Andel_1365066/2020
	* Netherlands/Helmond_1363548/2020

* MSHS Clinical Microbiology Laboratories
	* USA/NY1-PV08001/2020

* Massachusetts Department of Public Health
	* USA/MA1/2020

* Monash Medical Centre
	* Australia/VIC01/2020

* NHC Key laboratory of Enteric Pathogenic Microbiology, Institute of Pathogenic Microbiology
	* Jiangsu/JS01/2020
	* Jiangsu/JS02/2020
	* Jiangsu/JS03/2020

* National Centre for Infectious Diseases
	* Singapore/12/2020
	* Singapore/13/2020
	* Singapore/14/2020
	* Singapore/3/2020
	* Singapore/4/2020

* National Influenza Center - National Institute of Hygiene and Epidemiology (NIHE)
	* Vietnam/VR03-38142/2020

* National Influenza Centre, National Public Health Laboratory, Kathmandu, Nepal
	* Nepal/61/2020

* National Institute for Viral Disease Control and Prevention, China CDC
	* Beijing/IVDC-BJ-005/2020
	* Chongqing/IVDC-CQ-001/2020
	* Henan/IVDC-HeN-002/2020
	* Jiangsu/IVDC-JS-001/2020
	* Jiangxi/IVDC-JX-002/2020
	* Shandong/IVDC-SD-001/2020
	* Shanghai/IVDC-SH-001/2020
	* Sichuan/IVDC-SC-001/2020
	* Wuhan/IVDC-HB-01/2019
	* Wuhan/IVDC-HB-04/2020
	* Wuhan/IVDC-HB-05/2019
	* Yunnan/IVDC-YN-003/2020

* National Public Health Laboratory
	* Singapore/11/2020

* National Public Health Laboratory, National Centre for Infectious Diseases
	* Singapore/10/2020
	* Singapore/7/2020
	* Singapore/8/2020
	* Singapore/9/2020

* Pathology Queensland
	* Australia/QLD01/2020
	* Australia/QLD02/2020
	* Australia/QLD03/2020
	* Australia/QLD04/2020
	* Australia/QLD09/2020

* Providence Regional Medical Center
	* USA/WA1/2020

* Public Health Ontario Laboratory
	* Canada/ON-PHL2445/2020
	* Canada/ON-VIDO-01/2020

* RIVM
	* Netherlands/Delft_1363424/2020
	* Netherlands/Diemen_1363454/2020
	* Netherlands/Loon_op_zand_1363512/2020
	* Netherlands/Oss_1363500/2020
	* NetherlandsL/Houten_1363498/2020

* Respiratory Virus Unit, Microbiology Services Colindale, Public Health England
	* England/01/2020
	* England/02/2020
	* England/09c/2020
	* England/200641094/2020
	* England/200690245/2020
	* England/200690300/2020
	* England/200690306/2020
	* England/200690756/2020
	* England/200940527/2020
	* England/200960041/2020
	* England/200960515/2020
	* England/200981386/2020
	* England/200990002/2020
	* England/200990006/2020
	* England/200990660/2020
	* England/200990723/2020
	* England/200990724/2020
	* England/200990725/2020
	* England/200991076/2020
	* England/201000003/2020
	* England/201040081/2020
	* England/201040141/2020

* Seattle Flu Study
	* USA/WA-S2/2020
	* USA/WA-S3/2020

* Second Hospital of Anhui Medical University
	* Hefei/2/2020

* Serology, Virology and OTDS Laboratories (SAViD), NSW Health Pathology Randwick
	* Sydney/3/2020

* Servicio Microbiologia. Hospital Clinico Universitario. Valencia.
	* Spain/Valencia1/2020
	* Spain/Valencia2/2020

* Shenzhen Key Laboratory of Pathogen and Immunity, National Clinical Research Center for Infectious Disease, Shenzhen Third People's Hospital
	* Shenzhen/SZTH-002/2020
	* Shenzhen/SZTH-003/2020
	* Shenzhen/SZTH-004/2020

* Shenzhen Third People's Hospital
	* Shenzhen/SZTH-001/2020

* Singapore General Hospital
	* Singapore/1/2020
	* Singapore/2/2020

* Singapore General Hospital, Molecular Laboratory, Division of Pathology
	* Singapore/5/2020
	* Singapore/6/2020

* Sorbonne Universite, Inserm et Assistance Publique-Hopitaux de Paris (Pitie Salpetriere)
	* France/IDF0626/2020

* South China Agricultural University
	* pangolin/Guandong/1/2019

* State Health Office Baden-Wuerttemberg
	* Germany/Baden-Wuerttemberg-1/2020

* Taiwan Centers for Disease Control
	* Taiwan/3/2020
	* Taiwan/4/2020

* Texas Department of State Health Services
	* USA/TX1/2020

* The Central Hospital Of Wuhan
	* Wuhan/HBCDC-HB-02/2020

* The National Institute of Public Health Center for Epidemiology and Microbiology
	* CzechRepublic/951/2020

* The University of Hong Kong - Shenzhen Hospital
	* Shenzhen/HKU-SZ-002/2020
	* Shenzhen/HKU-SZ-005/2020

* Tianmen Center for Disease Control and Prevention
	* Tianmen/HBCDC-HB-07/2020

* UCD National Virus Reference Laboratory
	* Ireland/COR-20134/2020

* UW Virology Lab
	* USA/WA-UW15/2020
	* USA/WA-UW16/2020
	* USA/WA-UW17/2020
	* USA/WA-UW18/2020
	* USA/WA-UW19/2020
	* USA/WA-UW20/2020
	* USA/WA-UW21/2020
	* USA/WA11-UW7/2020
	* USA/WA12-UW8/2020
	* USA/WA13-UW9/2020
	* USA/WA14-UW10/2020
	* USA/WA15-UW11/2020
	* USA/WA16-UW12/2020
	* USA/WA17-UW13/2020
	* USA/WA18-UW14/2020

* Union Hospital of Tongji Medical College, Huazhong University of Science and Technology
	* Wuhan/HBCDC-HB-03/2020
	* Wuhan/HBCDC-HB-04/2020

* Unknown
	* Netherlands/Coevorden_1363618/2020

* Valley Medical Center
	* USA/WA8-UW5/2020

* Virology Department, Sheffield Teaching Hospitals NHS Foundation Trust
	* England/Sheff01/2020
	* England/Sheff02/2020

* Virology Unit, Institut Pasteur du Cambodge.
	* Cambodia/0012/2020

* WA State Department of Health
	* USA/WA1-A12/2020

* Wales Specialist Virology Centre
	* Wales/PHW03/2020
	* Wales/PHW05/2020
	* Wales/PHW1/2020
	* Wales/PHW2/2020

* Washington State Department of Health
	* USA/WA1-F6/2020
	* USA/WA2/2020

* Washington State Public Health Lab
	* USA/WA4-UW2/2020
	* USA/WA6-UW3/2020
	* USA/WA7-UW4/2020

* Weifang Center for Disease Control and Prevention
	* China/WF0001/2020
	* China/WF0002/2020
	* China/WF0003/2020
	* China/WF0004/2020
	* China/WF0006/2020
	* China/WF0009/2020
	* China/WF0012/2020
	* China/WF0014/2020
	* China/WF0015/2020
	* China/WF0016/2020
	* China/WF0017/2020
	* China/WF0018/2020
	* China/WF0019/2020
	* China/WF0020/2020
	* China/WF0021/2020
	* China/WF0023/2020
	* China/WF0024/2020
	* China/WF0026/2020
	* China/WF0028/2020
	* China/WF0029/2020

* West of Scotland Specialist Virology Centre, NHSGGC
	* Scotland/CVR01/2020
	* Scotland/CVR02/2020
	* Scotland/CVR03/2020
	* Scotland/CVR04/2020
	* Scotland/CVR05/2020

* Wisconsin Department of Health Services
	* USA/WI1/2020

* Wuhan Fourth Hospital
	* Wuhan/WH05/2020

* Wuhan Institute of Virology, Chinese Academy of Sciences
	* bat/Yunnan/RaTG13/2013

* Wuhan Jinyintan Hospital
	* Wuhan/HBCDC-HB-01/2019
	* Wuhan/HBCDC-HB-02/2019
	* Wuhan/HBCDC-HB-03/2019
	* Wuhan/HBCDC-HB-04/2019
	* Wuhan/WIV02/2019
	* Wuhan/WIV04/2019
	* Wuhan/WIV05/2019
	* Wuhan/WIV06/2019
	* Wuhan/WIV07/2019

* Wuhan Lung Hospital
	* Wuhan/HBCDC-HB-06/2020

* Yongchuan District Center for Disease Control and Prevention
	* Chongqing/YC01/2020

* Zhejiang Provincial Center for Disease Control and Prevention
	* Zhejiang/WZ-01/2020
	* Zhejiang/WZ-02/2020

* Zhongxian Center for Disease Control and Prevention
	* Chongqing/ZX01/2020

* Andersen Lab, The Scripps Research Institute
	* USA/CA-PC101P/2020

* Arizona Department of Health Services
	* USA/AZ1/2020

* Auckland Hospital
	* NewZealand/01/2020

* BCCDC Public Health Laboratory
	* Canada/BC_02421/2020
	* Canada/BC_13297/2020
	* Canada/BC_17397/2020
	* Canada/BC_25211/2020
	* Canada/BC_35720/2020
	* Canada/BC_37_0-2/2020
	* Canada/BC_40860/2020
	* Canada/BC_41851/2020
	* Canada/BC_64686/2020
	* Canada/BC_65034/2020
	* Canada/BC_66353/2020
	* Canada/BC_69243/2020
	* Canada/BC_78548/2020
	* Canada/BC_83109/2020
	* Canada/BC_83163/2020

* Bamrasnaradura Hospital
	* Nonthaburi/61/2020
	* Nonthaburi/74/2020

* Beijing Institute of Microbiology and Epidemiology
	* pangolin/Guangdong/P2S/2019
	* pangolin/Guangxi/P1E/2017
	* pangolin/Guangxi/P2V/2017
	* pangolin/Guangxi/P3B/2017
	* pangolin/Guangxi/P4L/2017
	* pangolin/Guangxi/P5E/2017
	* pangolin/Guangxi/P5L/2017

* Bundeswehr Institute of Microbiology
	* Germany/BavPat2/2020
	* Germany/BavPat3/2020

* CNR Virus des Infections Respiratoires - France SUD
	* France/RA739/2020

* CR&WISCO GENERAL HOSPITAL
	* Wuhan/HBCDC-HB-05/2020

* California Department of Health
	* USA/CA3/2020
	* USA/CA4/2020
	* USA/CA5/2020

* California Department of Public Health
	* USA/CA-CDPH-UC1/2020
	* USA/CA-CDPH-UC2/2020
	* USA/CA-CDPH-UC3/2020
	* USA/CA-CDPH-UC4/2020
	* USA/CA-CDPH-UC5/2020
	* USA/CA-CDPH-UC6/2020
	* USA/CA-CDPH-UC7/2020
	* USA/CA-CDPH-UC8/2020
	* USA/CA-CDPH-UC9/2020
	* USA/CA1/2020
	* USA/CA2/2020
	* USA/CA6/2020
	* USA/CA7/2020
	* USA/CA8/2020
	* USA/CA9/2020
	* USA/UC-CDPH-UC11/2020

* Center of Medical Microbiology, Virology, and Hospital Hygiene, University of Duesseldorf
	* Germany/NRW-01/2020
	* Germany/NRW-011/2020
	* Germany/NRW-02-1/2020
	* Germany/NRW-03/2020
	* Germany/NRW-04/2020
	* Germany/NRW-05/2020
	* Germany/NRW-06/2020
	* Germany/NRW-07/2020
	* Germany/NRW-08/2020
	* Germany/NRW-09/2020
	* Germany/NRW-10/2020

* Centers for Disease Control, R.O.C. (Taiwan)
	* Taiwan/2/2020

* Centre Hositalier Universitaire de Rouen Laboratoire de Virologie
	* France/N1620/2020

* Centre Hospitalier Compiegne Laboratoire de Biologie
	* France/HF1795/2020
	* France/HF1805/2020
	* France/HF1870/2020
	* France/HF1871/2020
	* France/HF1986/2020
	* France/HF1988/2020
	* France/HF1989/2020
	* France/HF1993/2020
	* France/HF1995/2020
	* France/HF2151/2020
	* France/HF2174/2020

* Centre Hospitalier Regional Universitaire de Nantes Laboratoire de Virologie
	* France/PL1643/2020

* Centre Hospitalier Rene Dubois Laboratoire de Microbiologie - Bat A
	* France/IDF1980/2020

* Centre for Human and Zoonotic Virology (CHAZVY), College of Medicine University of Lagos/Lagos University Teaching Hospital (LUTH), part of the Laboratory Network of the Nigeria Centre for Disease Control (NCDC)
	* Nigeria/Lagos01/2020

* Centre for Infectious Diseases and Microbiology - Public Health
	* Australia/NSW10/2020
	* Australia/NSW12/2020
	* Australia/NSW13/2020
	* Australia/NSW14/2020

* Centre for Infectious Diseases and Microbiology Laboratory Services
	* Australia/NSW01/2020
	* Australia/NSW02/2020
	* Australia/NSW05/2020
	* Australia/NSW06/2020
	* Australia/NSW07/2020
	* Australia/NSW08/2020
	* Australia/NSW09/2020

* Centre for Infectious Diseases and Microbiology- Public Health
	* Australia/NSW11/2020

* Centro Hospital do Porto, E.P.E. - H. Geral de Santo Antonio
	* Portugal/CV62/2020

* Centro Hospitalar e Universitario de Sao Joao, Porto
	* Portugal/CV63/2020

* Charite Universitatsmedizin Berlin, Institute of Virology; Institut fur Mikrobiologie der Bundeswehr, Munich
	* Germany/BavPat1/2020

* Clinica Alemana de Santiago, Chile
	* Chile/Santiago-1/2020

* Clinica Santa Maria, Santiago, Chile
	* Chile/Santiago-2/2020

* Department of Infectious Diseases, Istituto Superiore di Sanita, Roma , Italy
	* Italy/CDG1/2020

* Department of Infectious Diseases, Istituto Superiore di Sanita, Rome, Italy
	* Italy/SPL1/2020

* Department of Infectious and Tropical Diseases, Bichat Claude Bernard Hospital, Paris
	* France/IDF0372-isl/2020
	* France/IDF0372/2020
	* France/IDF0373/2020
	* France/IDF0386-islP1/2020
	* France/IDF0386-islP3/2020
	* France/IDF0515-isl/2020
	* France/IDF0515/2020
	* France/IDF0571/2020

* Department of Internal Medicine, Triemli Hospital
	* Switzerland/1000477102/2020
	* Switzerland/1000477377/2020

* Department of Laboratory Medicine, National Taiwan University Hospital
	* Taiwan/NTU01/2020
	* Taiwan/NTU02/2020
	* Taiwan/NTU03/2020

* Department of Microbiology, Institute for Viral Diseases, College of Medicine, Korea University
	* SouthKorea/KUMC01/2020
	* SouthKorea/KUMC02/2020
	* SouthKorea/KUMC04/2020
	* SouthKorea/KUMC06/2020

* Department of Pathology, Toshima Hospital
	* Japan/TK/20-31-3/2020

* Department of Virology III, National Institute of Infectious Diseases
	* Japan/AI/I-004/2020

* Department of Virology and Immunology, University of Helsinki and Helsinki University Hospital, Huslab Finland
	* Finland/FIN-114/2020
	* Finland/FIN-266/2020
	* Finland/FIN-274/2020
	* Finland/FIN-313/2020
	* Finland/FIN-318/2020
	* Finland/FIN-455/2020
	* Finland/FIN-508/2020
	* Finland/FIN01032020/2020
	* Finland/FIN03032020A/2020
	* Finland/FIN03032020B/2020
	* Finland/FIN03032020C/2020

* Department of Virus and Microbiological Special diagnostics, Statens Serum Institut, Copenhagen, Denmark.
	* Denmark/SSI-101/2020
	* Denmark/SSI-104/2020

* Department of microbiology laboratory,Anhui Provincial Center for Disease Control and Prevention
	* Anhui/SZ005/2020

* Dept. of Pathology, National Institute of Infectious Diseases
	* Japan/NA-20-05-1/2020
	* Japan/OS-20-07-1/2020

* Dept. of Virology III, National Institute of Infectious Diseases
	* Japan/KY-V-029/2020
	* Japan/TY-WK-012/2020
	* Japan/TY-WK-501/2020
	* Japan/TY-WK-521/2020

* Dienst Gezondheid & Jeugd Zuid-Holland Zuid
	* Netherlands/Hardinxveld_Giessendam_1364806/2020

* Division of Infectious Diseases, Department of Internal Medicine, Korea University College of Medicine
	* SouthKorea/KUMC03/2020
	* SouthKorea/KUMC05/2020

* Division of Infectious Diseases, University Hospital Zurich
	* Switzerland/1000477796/2020
	* Switzerland/1000477797/2020
	* Switzerland/1000477806/2020

* Division of Viral Diseases, Center for Laboratory Control of Infectious Diseases, Korea Centers for Diseases Control and Prevention
	* SouthKorea/KCDC05/2020
	* SouthKorea/KCDC06/2020
	* SouthKorea/KCDC07/2020
	* SouthKorea/KCDC12/2020
	* SouthKorea/KCDC24/2020

* Dutch COVID-19 response team
	* Netherlands/Flevoland/1/2020
	* Netherlands/Gelderland/1/2020
	* Netherlands/Gelderland/2/2020
	* Netherlands/Gelderland/3/2020
	* Netherlands/Gelderland_1/2020
	* Netherlands/Limburg/7/2020
	* Netherlands/Limburg_2/2020
	* Netherlands/Limburg_3/2020
	* Netherlands/Limburg_4/2020
	* Netherlands/Limburg_5/2020
	* Netherlands/Limburg_6/2020
	* Netherlands/NA/1/2020
	* Netherlands/NA/10/2020
	* Netherlands/NA/11/2020
	* Netherlands/NA/12/2020
	* Netherlands/NA/13/2020
	* Netherlands/NA/14/2020
	* Netherlands/NA/15/2020
	* Netherlands/NA/16/2020
	* Netherlands/NA/17/2020
	* Netherlands/NA/18/2020
	* Netherlands/NA/19/2020
	* Netherlands/NA/2/2020
	* Netherlands/NA/20/2020
	* Netherlands/NA/21/2020
	* Netherlands/NA/22/2020
	* Netherlands/NA/23/2020
	* Netherlands/NA/24/2020
	* Netherlands/NA/25/2020
	* Netherlands/NA/26/2020
	* Netherlands/NA/27/2020
	* Netherlands/NA/28/2020
	* Netherlands/NA/29/2020
	* Netherlands/NA/30/2020
	* Netherlands/NA/31/2020
	* Netherlands/NA/32/2020
	* Netherlands/NA/33/2020
	* Netherlands/NA/34/2020
	* Netherlands/NA/35/2020
	* Netherlands/NA/4/2020
	* Netherlands/NA/5/2020
	* Netherlands/NA/6/2020
	* Netherlands/NA/7/2020
	* Netherlands/NA/8/2020
	* Netherlands/NA/9/2020
	* Netherlands/NoordBrabant/41/2020
	* Netherlands/NoordBrabant/42/2020
	* Netherlands/NoordBrabant/44/2020
	* Netherlands/NoordBrabant/45/2020
	* Netherlands/NoordBrabant/46/2020
	* Netherlands/NoordBrabant/47/2020
	* Netherlands/NoordBrabant/48/2020
	* Netherlands/NoordBrabant/49/2020
	* Netherlands/NoordBrabant/51/2020
	* Netherlands/NoordBrabant/52/2020
	* Netherlands/NoordBrabant/53/2020
	* Netherlands/NoordBrabant/54/2020
	* Netherlands/NoordBrabant/55/2020
	* Netherlands/NoordBrabant/56/2020
	* Netherlands/NoordBrabant/57/2020
	* Netherlands/NoordBrabant/58/2020
	* Netherlands/NoordBrabant/59/2020
	* Netherlands/NoordBrabant/60/2020
	* Netherlands/NoordBrabant/61/2020
	* Netherlands/NoordBrabant/62/2020
	* Netherlands/NoordBrabant/63/2020
	* Netherlands/NoordBrabant/64/2020
	* Netherlands/NoordBrabant/65/2020
	* Netherlands/NoordBrabant/66/2020
	* Netherlands/NoordBrabant/67/2020
	* Netherlands/NoordBrabant/68/2020
	* Netherlands/NoordBrabant_1/2020
	* Netherlands/NoordBrabant_10/2020
	* Netherlands/NoordBrabant_11/2020
	* Netherlands/NoordBrabant_12/2020
	* Netherlands/NoordBrabant_13/2020
	* Netherlands/NoordBrabant_14/2020
	* Netherlands/NoordBrabant_15/2020
	* Netherlands/NoordBrabant_16/2020
	* Netherlands/NoordBrabant_17/2020
	* Netherlands/NoordBrabant_18/2020
	* Netherlands/NoordBrabant_19/2020
	* Netherlands/NoordBrabant_2/2020
	* Netherlands/NoordBrabant_20/2020
	* Netherlands/NoordBrabant_21/2020
	* Netherlands/NoordBrabant_22/2020
	* Netherlands/NoordBrabant_23/2020
	* Netherlands/NoordBrabant_24/2020
	* Netherlands/NoordBrabant_25/2020
	* Netherlands/NoordBrabant_26/2020
	* Netherlands/NoordBrabant_27/2020
	* Netherlands/NoordBrabant_28/2020
	* Netherlands/NoordBrabant_29/2020
	* Netherlands/NoordBrabant_3/2020
	* Netherlands/NoordBrabant_30/2020
	* Netherlands/NoordBrabant_31/2020
	* Netherlands/NoordBrabant_32/2020
	* Netherlands/NoordBrabant_33/2020
	* Netherlands/NoordBrabant_34/2020
	* Netherlands/NoordBrabant_35/2020
	* Netherlands/NoordBrabant_36/2020
	* Netherlands/NoordBrabant_37/2020
	* Netherlands/NoordBrabant_38/2020
	* Netherlands/NoordBrabant_39/2020
	* Netherlands/NoordBrabant_4/2020
	* Netherlands/NoordBrabant_5/2020
	* Netherlands/NoordBrabant_6/2020
	* Netherlands/NoordHolland/3/2020
	* Netherlands/NoordHolland_1/2020
	* Netherlands/NoordHolland_2/2020
	* Netherlands/Overijssel_1/2020
	* Netherlands/Overijssel_2/2020
	* Netherlands/Utrecht/17/2020
	* Netherlands/Utrecht/18/2020
	* Netherlands/Utrecht/19/2020
	* Netherlands/Utrecht_1/2020
	* Netherlands/Utrecht_10/2020
	* Netherlands/Utrecht_11/2020
	* Netherlands/Utrecht_12/2020
	* Netherlands/Utrecht_13/2020
	* Netherlands/Utrecht_14/2020
	* Netherlands/Utrecht_15/2020
	* Netherlands/Utrecht_16/2020
	* Netherlands/Utrecht_2/2020
	* Netherlands/Utrecht_3/2020
	* Netherlands/Utrecht_4/2020
	* Netherlands/Utrecht_5/2020
	* Netherlands/Utrecht_6/2020
	* Netherlands/Utrecht_7/2020
	* Netherlands/Utrecht_8/2020
	* Netherlands/ZuidHolland/25/2020
	* Netherlands/ZuidHolland/26/2020
	* Netherlands/ZuidHolland/27/2020
	* Netherlands/ZuidHolland/28/2020
	* Netherlands/ZuidHolland/29/2020
	* Netherlands/ZuidHolland/30/2020
	* Netherlands/ZuidHolland/31/2020
	* Netherlands/ZuidHolland_1/2020
	* Netherlands/ZuidHolland_10/2020
	* Netherlands/ZuidHolland_11/2020
	* Netherlands/ZuidHolland_13/2020
	* Netherlands/ZuidHolland_14/2020
	* Netherlands/ZuidHolland_15/2020
	* Netherlands/ZuidHolland_16/2020
	* Netherlands/ZuidHolland_17/2020
	* Netherlands/ZuidHolland_18/2020
	* Netherlands/ZuidHolland_19/2020
	* Netherlands/ZuidHolland_2/2020
	* Netherlands/ZuidHolland_20/2020
	* Netherlands/ZuidHolland_21/2020
	* Netherlands/ZuidHolland_22/2020
	* Netherlands/ZuidHolland_23/2020
	* Netherlands/ZuidHolland_24/2020
	* Netherlands/ZuidHolland_5/2020
	* Netherlands/ZuidHolland_6/2020
	* Netherlands/ZuidHolland_7/2020
	* Netherlands/ZuidHolland_8/2020
	* Netherlands/ZuidHolland_9/2020

* ErasmusMC
	* Netherlands/Nieuwendijk_1363582/2020
	* Netherlands/Rotterdam_1363790/2020

* Foundation Elisabeth-Tweesteden Ziekenhuis
	* Netherlands/Tilburg_1363354/2020
	* Netherlands/Tilburg_1364286/2020

* Foundation Pamm
	* Netherlands/Berlicum_1363564/2020

* Fujian Center for Disease Control and Prevention
	* Fujian/13/2020
	* Fujian/8/2020

* General Hospital of Central Theater Command of People's Liberation Army of China
	* Wuhan/WH01/2019
	* Wuhan/WH02/2019
	* Wuhan/WH03/2020
	* Wuhan/WH04/2020

* Gorgas Memorial Institute for Health Studies
	* Panama/328677/2020

* Guangdong Provincial Center for Diseases Control and Prevention; Guangdong Provincial Public Health
	* Foshan/20SF207/2020
	* Foshan/20SF210/2020
	* Foshan/20SF211/2020
	* Guangdong/20SF012/2020
	* Guangdong/20SF013/2020
	* Guangdong/20SF014/2020
	* Guangdong/20SF025/2020
	* Guangdong/20SF028/2020
	* Guangdong/20SF040/2020

* Guangdong Provincial Center for Diseases Control and Prevention; Guangdong Provinical Public Health
	* Guangdong/20SF174/2020
	* Guangzhou/20SF206/2020

* Guangdong Provincial Center for Diseases Control and Prevention;Guangdong Provincial Institute of Public Health
	* Guangdong/20SF201/2020

* Guangdong Provincial Institution of Public Health, Guangdong Provinical Center for Disease Control and Prevention
	* Guangdong/2020XN4239-P0034/2020
	* Guangdong/2020XN4243-P0035/2020
	* Guangdong/2020XN4273-P0036/2020
	* Guangdong/2020XN4276-P0037/2020
	* Guangdong/2020XN4291-P0038/2020
	* Guangdong/2020XN4373-P0039/2020
	* Guangdong/2020XN4433-P0040/2020
	* Guangdong/2020XN4448-P0002/2020
	* Guangdong/2020XN4459-P0041/2020
	* Guangdong/2020XN4475-P0042/2020
	* Guangdong/DG-S2-P0054/2020
	* Guangdong/DG-S41-P0056/2020
	* Guangdong/DG-S6-P0055/2020
	* Guangdong/DG-S9-P0045/2020
	* Guangdong/FS-S29-P0051/2020
	* Guangdong/FS-S30-P0052/2020
	* Guangdong/FS-S34-P0015/2020
	* Guangdong/FS-S42-P0046/2020
	* Guangdong/FS-S48-P0047/2020
	* Guangdong/FS-S50-P0053/2020
	* Guangdong/GD2020012-P0022/2020
	* Guangdong/GD2020016-P0011/2020
	* Guangdong/GD2020080-P0010/2020
	* Guangdong/GD2020085-P0043/2020
	* Guangdong/GD2020086-P0021/2020
	* Guangdong/GD2020087-P0008/2020
	* Guangdong/GD2020115-P0009/2020
	* Guangdong/GD2020134-P0031/2020
	* Guangdong/GD2020139-P0007/2020
	* Guangdong/GD2020227-P0029/2020
	* Guangdong/GD2020233-P0027/2020
	* Guangdong/GD2020234-P0023/2020
	* Guangdong/GD2020241-P0013/2020
	* Guangdong/GD2020246-P0028/2020
	* Guangdong/GD2020258-P0018/2020
	* Guangdong/GDFS2020052-P0025/2020
	* Guangdong/GDFS2020054-P0005/2020
	* Guangdong/GDFS2020056-P0044/2020
	* Guangdong/GDFS2020127-P0026/2020
	* Guangdong/GDSZ202004-P0004/2020
	* Guangdong/GDSZ202008-P0020/2020
	* Guangdong/GDSZ202009-P0032/2020
	* Guangdong/GDSZ202013-P0014/2020
	* Guangdong/GDSZ202015-P0019/2020
	* Guangdong/GZ-S6-P0050/2020
	* Guangdong/JM-S1-P0062/2020
	* Guangdong/MM-S1-P0048/2020
	* Guangdong/SZ-N128-P0057/2020
	* Guangdong/SZ-N59-P0049/2020
	* Guangdong/ZH-N22-P0059/2020
	* Guangdong/ZH-S33-P0058/2020
	* Guangdong/ZQ-S2-P0061/2020
	* Guangdong/ZS-S6-P0060/2020

* HUS Diagnostiikkakeskus, Hallinto
	* Finland/FIN-25/2020

* Hangzhou Center for Disease Control and Prevention
	* Hangzhou/HZCDC0001/2020

* Hangzhou Center for Disease and Control Microbiology Lab
	* Hangzhou/HZ-1/2020

* Harborview Medical Center
	* USA/WA3-UW1/2020
	* USA/WA9-UW6/2020

* Hong Kong Department of Health
	* HongKong/VB20024950/2020
	* HongKong/VB20026565/2020
	* HongKong/VM20001061/2020
	* HongKong/case42_VM20002493/2020
	* HongKong/case48_VM20002507/2020
	* HongKong/case49_VM20002508/2020
	* HongKong/case52_VM20002582/2020
	* HongKong/case78_VM20002849/2020
	* HongKong/case85_VM20002868/2020
	* HongKong/case90_VM20002907/2020

* Hopital Instruction des Armees - BEGIN
	* France/IDF2075/2020

* Hopital Robert Debre Laboratoire de Virologie
	* France/GE1973/2020
	* France/GE1977/2020

* Hopitaux universitaires de Geneve Laboratoire de Virologie
	* Switzerland/AG7120/2020
	* Switzerland/BE2536/2020
	* Switzerland/BE6651/2020
	* Switzerland/BS0914/2020
	* Switzerland/GE0199/2020
	* Switzerland/GE062072020
	* Switzerland/GE1402/2020
	* Switzerland/GE1422/2020
	* Switzerland/GE4135/2020
	* Switzerland/GE4984/2020
	* Switzerland/GE6679/2020
	* Switzerland/GE8102/2020
	* Switzerland/GR2988/2020
	* Switzerland/GR3043/2020
	* Switzerland/SZ1417/2020
	* Switzerland/TI2045/2020
	* Switzerland/VD0503/2020

* Hospital Israelita Albert Einstein
	* Brazil/SPBR-01/2020
	* Brazil/SPBR-02/2020
	* Brazil/SPBR-03/2020

* Hospital Sao Joaquim Beneficencia Portuguesa
	* Brazil/SPBR-04/2020
	* Brazil/SPBR-05/2020
	* Brazil/SPBR-06/2020

* Hospital de Talca, Chile
	* Chile/Talca-1/2020
	* Chile/Talca-2/2020

* IL Department of Public Health Chicago Laboratory
	* USA/IL1/2020
	* USA/IL2/2020

* INMI Lazzaro Spallanzani IRCCS
	* Italy/INMI1-cs/2020
	* Italy/INMI1-isl/2020

* Indian Council of Medical Research - National Institute of Virology
	* India/1-27/2020

* Indian Council of Medical Research-National Institute of Virology
	* India/1-31/2020

* Institute of Pathogen Biology, Chinese Academy of Medical Sciences & Peking Union Medical College
	* Wuhan/IPBCAMS-WH-01/2019
	* Wuhan/IPBCAMS-WH-02/2019
	* Wuhan/IPBCAMS-WH-03/2019
	* Wuhan/IPBCAMS-WH-04/2019
	* Wuhan/IPBCAMS-WH-05/2020

* Institute of Viral Disease Control and Prevention, China CDC
	* Wuhan/IVDC-HB-envF13-20/2020
	* Wuhan/IVDC-HB-envF13-21/2020
	* Wuhan/IVDC-HB-envF13/2020
	* Wuhan/IVDC-HB-envF54/2020

* Instituto Nacional de Enfermedades Respiratorias
	* Mexico/CDMX-InDRE_01/2020

* Jingzhou Center for Disease Control and Prevention
	* Jingzhou/HBCDC-HB-01/2020

* KU Leuven, Clincal and Epidemiological Virology
	* Belgium/BM-03012/2020

* KU Leuven, Clinical and Epidemiological Virology
	* Belgium/BA-02291/2020
	* Belgium/BC-03016/2020
	* Belgium/GHB-03021/2020
	* Belgium/QKJ-03015/2020
	* Belgium/SH-03014/2020
	* Belgium/VAG-03013/2020
	* Belgium/VLM-03011/2020

* Klinik Hirslanden Zurich
	* Switzerland/1000477757/2020

* Korea Centers for Disease Control & Prevention (KCDC) Center for Laboratory Control of Infectious Diseases Division of Viral Diseases
	* SouthKorea/KCDC03/2020

* LACEN RJ - Laboratorio Central de Saude Publica Noel Nutels
	* Brazil-RJ/314/2020

* LACEN/ES - Laboratorio Central de Saude Publica do Espirito Santo
	* Brazil/ES-225/2020

* Laboratoire National de Sante
	* Luxembourg/Lux1/2020

* Laboratoire de Virologie Institut de Virologie - INSERM U 1109 Hopitaux Universitaires de Strasbourg
	* France/GE1583/2020

* Laboratoire de Virologie, HUG
	* Switzerland/AG0361/2020
	* Switzerland/BL0902/2020
	* Switzerland/GE3121/2020
	* Switzerland/GE3895/2020
	* Switzerland/GE5373/2020
	* Switzerland/GE9586/2020
	* Switzerland/TI9486/2020
	* Switzerland/VD5615/2020

* Laboratorio Central de Saude Publica Professor Goncalo Moniz  LACEN/BA
	* Brazil/BA-312/2020

* Laboratorio di Microbiologia e Virologia, Universita Vita-Salute San Raffaele, Milano
	* Italy/UniSR1/2020

* Laboratory Medicine
	* Taiwan/CGMH-CGU-01/2020
	* Taiwan/CGMH-CGU-03/2020

* Laboratory of Molecular Virology, Pontificia Universidad Catolica de Chile
	* Chile/Santiago_op2d1/2020
	* Chile/Santiago_op3d1/2020
	* Chile/Santiago_op4d1/2020

* Lapland Central Hospital
	* Finland/1/2020

* MHC Brabant Zuidoost
	* Netherlands/Eindhoven_1363782/2020

* MHC Drente
	* Netherlands/Dalen_1363624/2020

* MHC Flevoland
	* Netherlands/Zeewolde_1365080/2020

* MHC Gooi & Vechtstreek
	* Netherlands/Blaricum_1364780/2020
	* Netherlands/Naarden_1364774/2020

* MHC Haaglanden
	* Netherlands/Nootdorp_1364222/2020

* MHC Hart voor Brabant
	* Netherlands/Oisterwijk_1364072/2020

* MHC Kennemerland
	* Netherlands/Haarlem_1363688/2020

* MHC Rotterdam-Rijnmond
	* Netherlands/Rotterdam_1364040/2020

* MHC Utrecht
	* Netherlands/Utrecht_1363564/2020
	* Netherlands/Utrecht_1363628/2020
	* Netherlands/Utrecht_1364066/2020

* MHC West-Brabant
	* Netherlands/Andel_1365066/2020
	* Netherlands/Helmond_1363548/2020

* MSHS Clinical Microbiology Laboratories
	* USA/NY1-PV08001/2020
	* USA/NY2-PV08100/2020

* Massachusetts Department of Public Health
	* USA/MA1/2020

* Minnesota Department of Health, Public Health Laboratory
	* USA/MN1-MDH1/2020
	* USA/MN2-MDH2/2020
	* USA/MN3-MDH3/2020

* Monash Medical Centre
	* Australia/VIC01/2020

* NHC Key laboratory of Enteric Pathogenic Microbiology, Institute of Pathogenic Microbiology
	* Jiangsu/JS01/2020
	* Jiangsu/JS02/2020
	* Jiangsu/JS03/2020

* NYU Langone Health
	* USA/NY-NYUMC1/2020

* National Centre for Infectious Diseases
	* Singapore/12/2020
	* Singapore/13/2020
	* Singapore/14/2020
	* Singapore/3/2020
	* Singapore/4/2020

* National Influenza Center - National Institute of Hygiene and Epidemiology (NIHE)
	* Vietnam/VR03-38142/2020

* National Influenza Centre, National Public Health Laboratory, Kathmandu, Nepal
	* Nepal/61/2020

* National Institute for Viral Disease Control and Prevention, China CDC
	* Beijing/IVDC-BJ-005/2020
	* Chongqing/IVDC-CQ-001/2020
	* Henan/IVDC-HeN-002/2020
	* Jiangsu/IVDC-JS-001/2020
	* Jiangxi/IVDC-JX-002/2020
	* Shandong/IVDC-SD-001/2020
	* Shanghai/IVDC-SH-001/2020
	* Sichuan/IVDC-SC-001/2020
	* Wuhan/IVDC-HB-01/2019
	* Wuhan/IVDC-HB-04/2020
	* Wuhan/IVDC-HB-05/2019
	* Yunnan/IVDC-YN-003/2020

* National Public Health Laboratory
	* Singapore/11/2020

* National Public Health Laboratory, National Centre for Infectious Diseases
	* Singapore/10/2020
	* Singapore/7/2020
	* Singapore/8/2020
	* Singapore/9/2020

* Pathology Queensland
	* Australia/QLD01/2020
	* Australia/QLD02/2020
	* Australia/QLD03/2020
	* Australia/QLD04/2020
	* Australia/QLD09/2020

* Providence Regional Medical Center
	* USA/WA1/2020

* Public Health Ontario Laboratory
	* Canada/ON-PHL2445/2020
	* Canada/ON-VIDO-01/2020

* R. G. Lugar Center for Public Health Research,  National Center for Disease Control and Public Health (NCDC) of Georgia.
	* Georgia/Tb-468/2020
	* Georgia/Tb-477/2020
	* Georgia/Tb-54/2020
	* Georgia/Tb-82/2020

* RIVM
	* Netherlands/Delft_1363424/2020
	* Netherlands/Diemen_1363454/2020
	* Netherlands/Loon_op_zand_1363512/2020
	* Netherlands/Oss_1363500/2020
	* NetherlandsL/Houten_1363498/2020

* Regional Virus Laboratory, Belfast
	* NorthernIreland/HSCNI01/2020

* Respiratory Virus Unit, Microbiology Services Colindale, Public Health England
	* England/01/2020
	* England/02/2020
	* England/09c/2020
	* England/200641094/2020
	* England/200690245/2020
	* England/200690300/2020
	* England/200690306/2020
	* England/200690756/2020
	* England/200940527/2020
	* England/200960041/2020
	* England/200960515/2020
	* England/200981386/2020
	* England/200990002/2020
	* England/200990006/2020
	* England/20099038206/2020
	* England/200990660/2020
	* England/200990723/2020
	* England/200990724/2020
	* England/200990725/2020
	* England/20099079106/2020
	* England/20099107406/2020
	* England/200991076/2020
	* England/201000003/2020
	* England/20100001406/2020
	* England/20100004706/2020
	* England/20100004806/2020
	* England/20100005406/2020
	* England/20100022706/2020
	* England/20100023206/2020
	* England/20100024006/2020
	* England/20100077906/2020
	* England/20100121006/2020
	* England/20100121007/2020
	* England/20100122106/2020
	* England/20100122107/2020
	* England/20102000106/2020
	* England/20102000206/2020
	* England/20102000306/2020
	* England/20102000506/2020
	* England/20102000906/2020
	* England/20102068506/2020
	* England/201040081/2020
	* England/201040141/2020
	* England/20110003506/2020

* Seattle Flu Study
	* USA/WA-S2/2020
	* USA/WA-S3/2020

* Second Hospital of Anhui Medical University
	* Hefei/2/2020

* Serology, Virology and OTDS Laboratories (SAViD), NSW Health Pathology Randwick
	* Australia/NSW03/2020

* Servicio Microbiologia, Hospital Clinico Universitario, Valencia
	* Spain/Valencia3/2020

* Servicio Microbiologia. Hospital Clinico Universitario. Valencia.
	* Spain/Valencia1/2020
	* Spain/Valencia2/2020

* Shandong Provincial Center for Disease Control and Prevention
	* Shandong/LY001/2020
	* Shandong/LY002/2020
	* Shandong/LY003/2020
	* Shandong/LY004/2020
	* Shandong/LY005/2020
	* Shandong/LY006/2020
	* Shandong/LY007/2020
	* Shandong/LY008/2020

* Shenzhen Key Laboratory of Pathogen and Immunity, National Clinical Research Center for Infectious Disease, Shenzhen Third People's Hospital
	* Shenzhen/SZTH-002/2020
	* Shenzhen/SZTH-003/2020
	* Shenzhen/SZTH-004/2020

* Shenzhen Third People's Hospital
	* Shenzhen/SZTH-001/2020

* Singapore General Hospital
	* Singapore/1/2020
	* Singapore/2/2020

* Singapore General Hospital, Molecular Laboratory, Division of Pathology
	* Singapore/5/2020
	* Singapore/6/2020

* Sorbonne Universite, Inserm et Assistance Publique-Hopitaux de Paris (Pitie Salpetriere)
	* France/IDF0626/2020

* South China Agricultural University
	* pangolin/Guangdong/1/2019

* State Health Office Baden-Wuerttemberg
	* Germany/Baden-Wuerttemberg-1/2020

* State Key Laboratory for Diagnosis and Treatment of Infectious Diseases, National Clinical Research Center for Infectious Diseases, First Affiliated Hospital, Zhejiang University School of Medicine, Hangzhou, China. 310003
	* Hangzhou/ZJU-01/2020
	* Hangzhou/ZJU-05/2020

* State Key Laboratory of Respiratory Disease, National Clinical Research Center for Respiratory Disease, Guangzhou Institute of Respiratory Health, the First Affiliated Hospital of Guangzhou Medical University
	* Guangzhou/GZMU0014/2020
	* Guangzhou/GZMU0016/2020
	* Guangzhou/GZMU0030/2020
	* Guangzhou/GZMU0031/2020
	* Guangzhou/GZMU0042/2020
	* Guangzhou/GZMU0044/2020
	* Guangzhou/GZMU0047/2020
	* Guangzhou/GZMU0048/2020

* Tai Lung Veterinary Laboratory, Agriculture, Fisheries and Conservation Department
	* canine/HongKong/20-02756/2020

* Taiwan Centers for Disease Control
	* Taiwan/3/2020
	* Taiwan/4/2020

* Texas Department of State Health Services
	* USA/TX1/2020

* The Central Hospital Of Wuhan
	* Wuhan/HBCDC-HB-02/2020

* The National Institute of Public Health Center for Epidemiology and Microbiology
	* CzechRepublic/951/2020

* The University of Hong Kong - Shenzhen Hospital
	* Shenzhen/HKU-SZ-002/2020
	* Shenzhen/HKU-SZ-005/2020

* Tianmen Center for Disease Control and Prevention
	* Tianmen/HBCDC-HB-07/2020

* UCD National Virus Reference Laboratory
	* Ireland/COR-20134/2020
	* Ireland/Dublin-19072/2020
	* Ireland/Limerick-19933/2020
	* Ireland/Limerick-19934/2020
	* Ireland/Limerick-19935/2020

* UW Virology Lab
	* USA/WA-UW15/2020
	* USA/WA-UW16/2020
	* USA/WA-UW17/2020
	* USA/WA-UW18/2020
	* USA/WA-UW19/2020
	* USA/WA-UW20/2020
	* USA/WA-UW21/2020
	* USA/WA-UW22/2020
	* USA/WA-UW23/2020
	* USA/WA-UW24/2020
	* USA/WA-UW25/2020
	* USA/WA-UW26/2020
	* USA/WA-UW27/2020
	* USA/WA-UW28/2020
	* USA/WA-UW29/2020
	* USA/WA-UW30/2020
	* USA/WA-UW31/2020
	* USA/WA-UW32/2020
	* USA/WA-UW33/2020
	* USA/WA-UW34/2020
	* USA/WA-UW35/2020
	* USA/WA-UW40/2020
	* USA/WA-UW41/2020
	* USA/WA-UW42/2020
	* USA/WA-UW43/2020
	* USA/WA-UW44/2020
	* USA/WA-UW45/2020
	* USA/WA-UW46/2020
	* USA/WA-UW47/2020
	* USA/WA-UW48/2020
	* USA/WA-UW49/2020
	* USA/WA-UW50/2020
	* USA/WA-UW51/2020
	* USA/WA-UW52/2020
	* USA/WA-UW53/2020
	* USA/WA-UW54/2020
	* USA/WA-UW55/2020
	* USA/WA-UW56/2020
	* USA/WA-UW57/2020
	* USA/WA-UW58/2020
	* USA/WA-UW59/2020
	* USA/WA-UW60/2020
	* USA/WA-UW61/2020
	* USA/WA-UW62/2020
	* USA/WA-UW63/2020
	* USA/WA-UW64/2020
	* USA/WA-UW65/2020
	* USA/WA-UW66/2020
	* USA/WA-UW67/2020
	* USA/WA-UW68/2020
	* USA/WA-UW69/2020
	* USA/WA-UW70/2020
	* USA/WA-UW71/2020
	* USA/WA-UW72/2020
	* USA/WA-UW73/2020
	* USA/WA-UW74/2020
	* USA/WA-UW75/2020
	* USA/WA-UW76/2020
	* USA/WA11-UW7/2020
	* USA/WA12-UW8/2020
	* USA/WA13-UW9/2020
	* USA/WA14-UW10/2020
	* USA/WA15-UW11/2020
	* USA/WA16-UW12/2020
	* USA/WA17-UW13/2020
	* USA/WA18-UW14/2020

* Union Hospital of Tongji Medical College, Huazhong University of Science and Technology
	* Wuhan/HBCDC-HB-03/2020
	* Wuhan/HBCDC-HB-04/2020

* Unknown
	* France/BFC2094/2020
	* Netherlands/Coevorden_1363618/2020

* Utah Public Health Laboratory
	* USA/UPHL-01/2020
	* USA/UPHL-02/2020
	* USA/UPHL-03/2020
	* USA/UPHL-04/2020
	* USA/UPHL-05/2020
	* USA/UPHL-06/2020

* Valley Medical Center
	* USA/WA8-UW5/2020

* Viral Respiratory Lab, National Institute for Biomedical Research (INRB)
	* Congo/KN-13/2020

* Virology Department, Royal Infirmary of Edinburgh, NHS Lothian
	* Scotland/EDB003/2020
	* Scotland/EDB004/2020
	* Scotland/EDB005/2020
	* Scotland/EDB006/2020
	* Scotland/EDB007/2020
	* Scotland/EDB008/2020
	* Scotland/EDB009/2020
	* Scotland/EDB011/2020
	* Scotland/EDB012/2020
	* Scotland/EDB013/2020

* Virology Department, Sheffield Teaching Hospitals NHS Foundation Trust
	* England/Sheff01/2020
	* England/Sheff02/2020

* Virology Unit, Institut Pasteur du Cambodge.
	* Cambodia/0012/2020

* WA State Department of Health
	* USA/WA1-A12/2020

* WHO National Influenza Centre Russian Federation
	* Russia/StPetersburg-3524/2020

* Wales Specialist Virology Centre
	* Wales/PHW03/2020
	* Wales/PHW04/2020
	* Wales/PHW05/2020
	* Wales/PHW06/2020
	* Wales/PHW1/2020
	* Wales/PHW10/2020
	* Wales/PHW12/2020
	* Wales/PHW13/2020
	* Wales/PHW2/2020
	* Wales/PHW27/2020
	* Wales/PHW28/2020
	* Wales/PHW33/2020
	* Wales/PHW37/2020
	* Wales/PHW38/2020

* Washington State Department of Health
	* USA/WA1-F6/2020
	* USA/WA2/2020

* Washington State Public Health Lab
	* USA/WA4-UW2/2020
	* USA/WA6-UW3/2020
	* USA/WA7-UW4/2020

* Weifang Center for Disease Control and Prevention
	* China/WF0001/2020
	* China/WF0002/2020
	* China/WF0003/2020
	* China/WF0004/2020
	* China/WF0006/2020
	* China/WF0009/2020
	* China/WF0012/2020
	* China/WF0014/2020
	* China/WF0015/2020
	* China/WF0016/2020
	* China/WF0017/2020
	* China/WF0018/2020
	* China/WF0019/2020
	* China/WF0020/2020
	* China/WF0021/2020
	* China/WF0023/2020
	* China/WF0024/2020
	* China/WF0026/2020
	* China/WF0028/2020
	* China/WF0029/2020

* West of Scotland Specialist Virology Centre, NHSGGC
	* Scotland/CVR01/2020
	* Scotland/CVR02/2020
	* Scotland/CVR03/2020
	* Scotland/CVR04/2020
	* Scotland/CVR05/2020
	* Scotland/CVR06/2020
	* Scotland/CVR07/2020
	* Scotland/CVR10/2020

* Wisconsin Department of Health Services
	* USA/WI1/2020

* Wuhan Fourth Hospital
	* Wuhan/WH05/2020

* Wuhan Institute of Virology, Chinese Academy of Sciences
	* bat/Yunnan/RaTG13/2013

* Wuhan Jinyintan Hospital
	* Wuhan/HBCDC-HB-01/2019
	* Wuhan/HBCDC-HB-02/2019
	* Wuhan/HBCDC-HB-03/2019
	* Wuhan/HBCDC-HB-04/2019
	* Wuhan/WIV02/2019
	* Wuhan/WIV04/2019
	* Wuhan/WIV05/2019
	* Wuhan/WIV06/2019
	* Wuhan/WIV07/2019

* Wuhan Lung Hospital
	* Wuhan/HBCDC-HB-06/2020

* Yongchuan District Center for Disease Control and Prevention
	* Chongqing/YC01/2020

* Zhejiang Provincial Center for Disease Control and Prevention
	* Zhejiang/WZ-01/2020
	* Zhejiang/WZ-02/2020

* Zhongxian Center for Disease Control and Prevention
	* Chongqing/ZX01/2020


```
