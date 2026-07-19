#!/bin/bash
cd /home/yhx/Code/health_physics/project/cell-mk/microtrack

echo "============================================"
echo "Am-241 siteRadius 扫描: 0.20 → 0.40 µm (21点, 单线程)"
echo "开始: $(date)"
echo "============================================"

for i in $(seq 0 20); do
    rd_nm=$((200 + i * 10))
    rd_um=$(awk "BEGIN{printf \"%.2f\", $rd_nm / 1000}")
    outfile="data/Am241/am241_phy_decay_rd${rd_nm}nm.root"

    if [ -f "$outfile" ] && [ $(stat -c%s "$outfile" 2>/dev/null || echo 0) -gt 100000 ]; then
        echo "[$((i+1))/21] rd=${rd_um}µm — 已存在, 跳过"
        continue
    fi

    echo "[$((i+1))/21] rd=${rd_um}µm (${rd_nm}nm) — 开始 $(date +%H:%M:%S)"

    sed -e "s|/mygeom/siteRadius.*|/mygeom/siteRadius  ${rd_um} um|" \
        -e "s|/run/numberOfThreads.*|/run/numberOfThreads 1|" \
        macro/run_am241_decay.mac > /tmp/am241_rd${rd_nm}.mac

    timeout 3600 ./build/microtrack /tmp/am241_rd${rd_nm}.mac > /tmp/am241_rd${rd_nm}.out 2>&1
    rc=$?

    if [ $rc -eq 0 ] && [ -f "$outfile" ] && [ $(stat -c%s "$outfile") -gt 100000 ]; then
        echo "  ✅ $(date +%H:%M:%S) ${outfile} ($(stat -c%s "$outfile") B)"
    else
        echo "  ❌ exit=$rc"
        rm -f "$outfile"
    fi
done

echo ""
echo "完成: $(date)"
echo "文件数: $(ls data/Am241/am241_phy_decay_rd*.root 2>/dev/null | wc -l)"
