#!/bin/bash

# 複数のTSPファイルを順番に実行するスクリプト
# 使用方法: ./run_multiple_tsp.sh

# 実行するTSPファイルのリスト（ディレクトリ名:ファイル名の形式）
# 形式: "ディレクトリ名:ファイル名"
TSP_FILES=(
    # "tsplib:rat575.tsp"
    # "tsplib:u724.tsp"
    "tsplib:vm1084.tsp"
    "tsplib:pcb1173.tsp"
    "tsplib:d1291.tsp"
    "tsplib:u1432.tsp"
    "tsplib:d1655.tsp"
    "tsplib:vm1748.tsp"
    "tsplib:u1817.tsp"
    "tsplib:u2152.tsp"
    "tsplib:pr2392.tsp"
    "tsplib:pcb3038.tsp"
    "tsplib:fnl4461.tsp"
    "vlsi:pla7397.tsp"
    "vlsi:rl5915.tsp"
    "vlsi:rbx711.tsp"
    "vlsi:xit1083.tsp"
    "vlsi:icw1483.tsp"
    "vlsi:djc1785.tsp"
    "vlsi:dcb2086.tsp"
    "vlsi:xpr2308.tsp"
    "vlsi:mlt2597.tsp"
    "vlsi:lsm2854.tsp"
    "vlsi:pia3056.tsp"
    "vlsi:fdp3256.tsp"
    "vlsi:ltb3729.tsp"
    "vlsi:bgb4355.tsp"
    "vlsi:xqd4966.tsp"
    "vlsi:fea5557.tsp"
    # ここに追加のファイルを追加できます
    # 形式: "ディレクトリ名:ファイル名"
)

# 実行パラメータ
POPULATION_SIZE=200
GENERATIONAL_MODEL="ER"
CHILDREN=100
SELECTION="greedy"
EAX_TYPE="EAX_5_AB"
TRIALS=30

# 結果ファイル名（各実行で追記される）
OUTPUT_FILE="result.md"

# 各TSPファイルに対して実行
for tsp_entry in "${TSP_FILES[@]}"; do
    # ディレクトリ名とファイル名を分割
    IFS=':' read -r directory tsp_file <<< "$tsp_entry"
    
    echo "=========================================="
    echo "Running: $directory/$tsp_file"
    echo "=========================================="
    
    # 指定されたディレクトリ内のファイルを指定
    make run/eax_tabu ARGS="--file $directory/$tsp_file --ps $POPULATION_SIZE --generational-model $GENERATIONAL_MODEL --children $CHILDREN --selection $SELECTION --eax-type $EAX_TYPE --trials $TRIALS --output $OUTPUT_FILE"
    
    # エラーチェック
    if [ $? -ne 0 ]; then
        echo "Error: Failed to run $directory/$tsp_file"
        echo "Continuing with next file..."
    fi
    
    echo ""
    echo "Completed: $directory/$tsp_file"
    echo ""
done

echo "=========================================="
echo "All TSP files processed!"
echo "Results saved to: $OUTPUT_FILE"
echo "=========================================="

