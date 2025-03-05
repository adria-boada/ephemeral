
# ephemeral

Cal instal·lar el mòdul `moments-popgen`.

```sh
# Crea un directori per a «virtual environments».
mkdir -p ~/.venvs_python/
cd ~/.venvs_python/
# Crea'l.
python3 -m venv moments-popgen
# Activa'l i instal·la els paquets necessaris.
source ~/.venvs_python/moments-popgen/bin/activate
pip install moments-popgen matplotlib

# Recorda activar l'ambient virtual "moments-popgen" abans d'utilitzar aquest paquet.

# Per utilitzar un dels ShellScripts, també cal la següent variable...
export MOMENTS_VIRTUAL_ENV_ACTIVATE="$HOME/.venvs_python/moments-popgen/bin/activate"
```

Per poder cridar al mòdul des de qualsevol «current working directory», cal ajustar
l'«environment variable» `PYTHONPATH`.

```sh
# Folder where 'ephemeral' is cloned. Could be copied to bashrc.
export PYTHONPATH="$PYTHONPATH:/home/user/folder"
```

