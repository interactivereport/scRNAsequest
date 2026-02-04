#!/usr/bin/env bash

if [[ -n "$1" ]]; then
    src="$(dirname $0)/src"
    if [[ ! -f "$src/.env" ]]; then
      echo "scPub: Please run ./install to setup necessary env variables"
      exit
    fi
    source $src/.env;eval $condaEnv
    if [[ ! -f "$HOME/.gitconfig" ]] || [[ $(grep "$(dirname $0)/.git" ~/.gitconfig | wc -l) == 0 ]];then # 
      git config --global --add safe.directory $(dirname $0)/.git
    fi

    env -i src="$src" a=$1 b=$2 c=$3 bash -c 'source $src/.env;eval $condaEnv;Rscript $src/publishCellDepot.R $src/ $a $b $c'

else
    echo "=============== scPub path/to/config/file"
    exePath=$(dirname $0)
    echo "'scPub' can be used to publish a scRNAseq result h5ad file to CellDepot."
    echo "Please check the example config files here: $exePath/demo/Pub2CellDepot_config.yml"
    echo "==============="
    exit 0
fi
