
# ephemeral

Crear _«virtual environment»_ amb els mòduls necessaris.

```sh
# Crea un directori per a «virtual environments».
mkdir -p ~/.venvs_python/
cd ~/.venvs_python/
python3 -m venv moments-popgen
# Activa el «v. env» i instal·la els paquets necessaris.
source ~/.venvs_python/moments-popgen/bin/activate
pip install moments-popgen matplotlib threadpoolctl

# Recorda activar l'ambient virtual "moments-popgen" abans d'utilitzar aquest paquet.

# Per utilitzar un dels ShellScripts, també cal la següent variable...
# (copiar a .bashrc per fer la variable permanent).
export MOMENTS_VIRTUAL_ENV_ACTIVATE="$HOME/.venvs_python/moments-popgen/bin/activate"
```

Per poder cridar al mòdul des de qualsevol _«current working directory»_, cal ajustar
l'_«environment variable»_ `PYTHONPATH`.

```sh
# Directori on 'ephemeral' s'ha clonat.
# (copiar a .bashrc per fer la variable permanent).
export PYTHONPATH="$PYTHONPATH:/home/user/folder"
```

Ara s'hauria de poder còrrer el paquet.

```sh
python3 -m ephemeral
# Les comandes estan ben documentades.
python3 -m ephemeral toSFS -h

# També es pot córrer en 'batches' a través d'un shellscript.
qsub -t $NUM_EXECUCIONS_INDEP ShellScripts/sge_moments_pipe.sh
```

