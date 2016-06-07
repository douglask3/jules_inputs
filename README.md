# Regridding to GC3 files

## Run JULES

- Compile jules. Set `fcm` path in `jules/compile_fules`   


```
    cd jules
    chmod +x compileJules
    ./compileJules ~/doukel/jules-vn4.5
```

- Run to

```
    cd jules
```
    - Edit namelists if required. Then:

```
    ./jules.exe
```

## Commands

```	
	./runNotebook.sh
```

Opens dir tree for all notebooks
```
	./runNotebook.sh regrid_clim.ipynb
```	

Open ipython notebook for regridding climate data


