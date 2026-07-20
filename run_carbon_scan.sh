#!/bin/bash
cd /home/yhx/Code/health_physics/project/cell-mk/microtrack

echo "============================================"
echo "C-12 siteRadius 扫描: 0.20 → 0.40 µm"
echo "21点, 2线程, 1000事件/点, timeout=3600s"
echo "开始: $(date)"
echo "============================================"

for i in $(seq 0 20); do
    rd_nm=$((200 + i * 10))
    rd_um=$(awk "BEGIN{printf \"%.2f\", $rd_nm / 1000}")
    outfile="data/Carbon/carbon_rd${rd_nm}nm.root"

    if [ -f "$outfile" ] && [ $(stat -c%s "$outfile" 2>/dev/null || echo 0) -gt 100000 ]; then
        echo "[$((i+1))/21] rd=${rd_um}µm — 已存在, 跳过"
        continue
    fi

    echo "[$((i+1))/21] rd=${rd_um}µm (${rd_nm}nm) — 开始 $(date +%H:%M:%S)"

    sed -e "s|/mygeom/siteRadius.*|/mygeom/siteRadius  ${rd_um} um|" \
        -e "s|/run/numberOfThreads.*|/run/numberOfThreads 4|" \
        macro/run_carbon.mac > /tmp/carbon_rd${rd_nm}.mac

    timeout 3600 ./build/microtrack /tmp/carbon_rd${rd_nm}.mac > /tmp/carbon_rd${rd_nm}.out 2>&1
    rc=$?

    if [ $rc -eq 0 ] && [ -f "$outfile" ] && [ $(stat -c%s "$outfile") -gt 100000 ]; then
        echo "  ✅ $(date +%H:%M:%S) $(stat -c%s "$outfile") B"
    else
        echo "  ❌ exit=$rc 大小=$(stat -c%s "$outfile" 2>/dev/null || echo 0) B"
        rm -f "$outfile"
    fi
done

echo ""
echo "完成: $(date)"
echo "文件数: $(ls data/Carbon/carbon_rd*.root 2>/dev/null | wc -l)"
ls data/Carbon/carbon_rd*.root 2>/dev/null
