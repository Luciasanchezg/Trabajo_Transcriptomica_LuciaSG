#### Creación del ficheros .gct ####

# ./Trabajo_6_generacion_gct.sh expression_HPB_ALL.csv HPB_ALL

INPUT=$1
NOMBRE=$2

# 1. Adición de la cabecera
echo '#1.2' > 'table_'$NOMBRE'.gct'

# 2. Creación de un archivo temporal con el formato adecuado
tail -n +2 $INPUT | tr "," "\t" > 'tmp_1_'$NOMBRE'.csv'

# 3. Adición del número de genes presentes
echo -n $(wc -l 'tmp_1_'$NOMBRE'.csv' | cut -d' ' -f 1)$'\t' >> 'table_'$NOMBRE'.gct'

# 4. Adición del número de muestras
echo '6' >> 'table_'$NOMBRE'.gct'

# 5. Adición de los títulos de las columnas
echo -n 'NAME'$'\t''Description' >> 'table_'$NOMBRE'.gct'
head -n 1 $INPUT | sed 's/,/\t/g' >> 'table_'$NOMBRE'.gct'

# 6. Ejecución del awk que introduce la columna NA
awk '{FS = "\t"; OFS = "\t"} { print $1, "NA" , $2, $3, $4, $5, $6, $7 }' > 'tmp_2_'$NOMBRE'.csv' 'tmp_1_'$NOMBRE'.csv'
cat 'tmp_2_'$NOMBRE'.csv' >> 'table_'$NOMBRE'.gct' 

# Eliminación de archivos temporales
rm 'tmp_1_'$NOMBRE'.csv'
rm 'tmp_2_'$NOMBRE'.csv'

