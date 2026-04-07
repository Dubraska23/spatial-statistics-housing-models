# 🏠 Spatial Statistics: Modelado Predictivo de Precios de Vivienda
### Análisis de Dependencia Espacial y Modelos SAR en R

Este proyecto aplica técnicas avanzadas de estadística espacial para entender y predecir el comportamiento de los precios inmobiliarios. A diferencia de los modelos lineales tradicionales, este enfoque captura la influencia de la ubicación y la proximidad (efecto vecindario) mediante econometría espacial.

---

## 🛡️ Visión de Gobierno de Datos y Calidad Espacial
El análisis de datos geográficos exige un control de calidad superior. En este proyecto se aplican los siguientes pilares de **Data Governance**:

1. **Integridad de Metadatos Geográficos:** Validación y estandarización del Sistema de Referencia de Coordenadas (CRS 4326 - WGS84) para asegurar que la geolocalización de los activos sea precisa y consistente.
2. **Control de Calidad en la Ingesta:** Protocolo de carga modular mediante el uso de carpetas estructuradas (`0_Data`), asegurando la trazabilidad de las funciones auxiliares y los microdatos originales.
3. **Validación de Dependencia (Moran’s I):** Implementación de pruebas de diagnóstico para detectar autocorrelación espacial en los residuos, garantizando que el modelo no presente sesgos por omisión de variables territoriales.
4. **Gobernanza de Modelos:** Comparativa de rendimiento entre Modelos Lineales (LM) y Modelos de Autorregresión Espacial (SAR) para asegurar que la decisión de negocio se base en el modelo estadísticamente más robusto.

---

## 🎯 Capacidades Técnicas
- **Análisis Exploratorio Espacial (ESDA):** Visualización de puntos y detección de clusters de precios.
- **Modelado SAR (Spatial Autoregressive):** Estimación del coeficiente Rho (ρ) para cuantificar la propagación del precio entre vecinos.
- **Capacidad Predictiva:** Evaluación de modelos mediante métricas de error (RMSE y R²) comparando datos de Train y Test (70/30).
- **Matrices de Pesos Espaciales:** Construcción de estructuras de vecindad para definir la interacción entre unidades geográficas.

---

## 🛠️ Stack Tecnológico
- **Lenguaje:** R
- **Librerías Espaciales:** `sf` (Simple Features), `spdep`, `spatialreg`.
- **Manipulación de Datos:** `tidyverse`, `dplyr`.
- **Reporting:** `R Markdown`.

---

## 📁 Estructura del Repositorio
Para garantizar la reproducibilidad, el proyecto mantiene la siguiente jerarquía:
- `/0_Data`: Contiene la base de datos `bbdd_12.csv` y scripts de librerías.
- `analisis_espacial_vivienda.Rmd`: Archivo principal de análisis y generación de reporte.

---

## 🚀 Cómo ejecutarlo
1. Clona el repositorio respetando la estructura de carpetas.
2. Abre el archivo `.Rmd` en RStudio.
3. Asegúrate de tener instalados los paquetes: `install.packages(c("sf", "spdep", "tidyverse", "spatialreg"))`.
4. Ejecuta el "Knit" para generar el informe final en HTML.

---
**Sobre la autora:** Especialista en **Data Governance** y **Business Analytics**. 
