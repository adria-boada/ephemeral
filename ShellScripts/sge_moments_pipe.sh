#! /bin/bash
#$ -cwd
#$ -V
#$ -j yes
#$ -S /bin/bash

# S'ha d'utilitzar la variable '$TASK_ID' per posar el número '$SGE_TASK_ID' al
# fitxer de sortida de l'«standard stream».
#$ -o $JOB_NAME.stdstream.job$JOB_ID.tsk$TASK_ID.log
#$ -N moments

# Funció per parar el guió amb un error.
error() {
	echo "ERROR:SGE_TASK_""$SGE_TASK_ID"":""$1" >&2
	echo "Terminating..." >&2
	exit "$2"
}

# Funció que s'utilitzarà per enviar informació detallada, per tal de
# diagnosticar problemes.
debug() {
	[ $debug -eq 0 ] && echo "DEBUG:SGE_TASK_""$SGE_TASK_ID"":""$@" >&2
}

# Retorna instruccions d'ús.
retorna_instruccions() {
	printf "Runs 'moments' optimization pipeline through 'qsub' of Sun Grid Engine.\n"
	printf "A multidimensional SFS (formatted as 'moments.Spectrum.tofile()') must be provided.\n\n"
	printf "Usage: %s -m|--model MODEL -s|--sfs 2D-SFS -o OUTFOLDER [--debug]\n" "$0"
	printf "\t%s\n" "--model: model requested for optimization."
	printf "\t%s\n" "$(echo "--model:" | tr "[:print:]" " ") See available models by running the module 'ephemeral optim --fit h' or in the script 'predefined_models.py'."
	printf "\t%s\n" "--sfs: the output SFS obtained from 'ephemeral toSFS'."
	printf "\t%s\n" "-o: a prefix that will be combined with the 'model' string to create an output directory."
	printf "\t%s\n" "--debug: print additional, detailed information helpful for diagnosing problems."
	exit 2
}

# Paràmetres
# ----------

# Variables predeterminades.
debug=1  # desactiva debug.

# Gestiona arguments.
while true ; do
	case "$1" in
		# Model triat que s'intentarà encaixar.
		-m | --model )
		[ -z "$2" ] && retorna_instruccions ; pop_model="$2" ; shift 2 ;;
		# Arxiu d'entrada multidimensional SFS.
		-s | --sfs )
		[ -z "$2" ] && retorna_instruccions ; arxiu_sfs="$2" ; shift 2 ;;
		# Crear carpeta output.
		-o )
		[ -z "$2" ] && retorna_instruccions ; outfold="$2" ; shift 2 ;;
		# Retorna missatges 'debug' si es passa aquesta opció.
		--debug ) debug=0 ; shift 1 ;;
		# Si no hi ha res hem acabat de llegir la línia de comandes.
		"" ) break ;;
		# Qualsevol altra cosa, retorna instruccions.
		* ) retorna_instruccions ; break ;;
	esac
done

# Revisa que "SGE_TASK_ID" funcioni i sigui un número.
if [[ ! "$SGE_TASK_ID" =~ ^[0-9]+$ ]] ; then
	msg="The parameter 'SGE_TASK_ID' is not provided (see 'man qsub',"
	msg="$msg option '-t')."
	error "$msg" 2
# Revisa que l'arxiu SFS existeix.
elif [ -z "$arxiu_sfs" ] ; then
	error "Please provide '-s SFS'." 2
elif [ ! -f "$arxiu_sfs" ] ; then
	error "The provided file '$arxiu_sfs' does not exist." 2
elif [ -z "$outfold" ] ; then
	error "Please provide '-o OUT-FOLDER'." 2
elif [ -z "$pop_model" ] ; then
	error "Please provide '-m MODEL'." 2
fi

# Es necessita especificar (exportar) les següents variables estàtiques
# a un arxiu separat. Podria ser .bashrc o un arxiu tmp (i llavors se'n fa
# 'source' cada vegada que es vulgui còrrer aquest guió).
if [ ! -f "$MOMENTS_VIRTUAL_ENV_ACTIVATE" ] ; then
	msg="Please, export an environment variable with the location of"
	msg="$msg the 'activate' shellscript of 'moments' module. The var."
	msg="$msg should be called 'MOMENTS_VIRTUAL_ENV_ACTIVATE'."
	error "$msg" 2
fi
# Try to find if the module 'ephemeral' is reachable with 'import':
python3 -c 'import ephemeral' 2> /dev/null || {
	msg="Please, update the environment variable 'PYTHONPATH' so it includes"
	msg="$msg the path to 'ephemeral' folder and export the variable."
	error "$msg" 2
}

# Compute the rest of the parameters.
carpeta_sortida="$outfold.$pop_model"
# Crea una carpeta per aglutinar els resultats de diferents tasques
# (execucions) paral·leles.
mkdir -p "$carpeta_sortida" && cd "$carpeta_sortida"

debug "Número de tasca SGE (execució) és $SGE_TASK_ID"
debug "Model seleccionat: $pop_model"
debug "Les etiquetes de les poblacions seleccionades (haurien de ser dues) són:"
debug "Nom de la carpeta de sortida hauria de ser '$carpeta_sortida'"
debug "L'arxiu SFS d'entrada és '$arxiu_sfs'"



source "$MOMENTS_VIRTUAL_ENV_ACTIVATE"
debug "Llançant guió d'optimizació de Python..."

python3 -m ephemeral optim "$arxiu_sfs" \
	--threads 1 \
	--fit-models "$pop_model" \
	--rounds 5 --replicates 10 20 20 50 50 \
	--manual-exec-i "$SGE_TASK_ID" \
	--max-iters 10 20 30 50 100 --fold-algor 3 2 2 2 1

debug "S'ha acabat l'anàlisis '$carpeta_sortida' de l'execució SGE $SGE_TASK_ID"

