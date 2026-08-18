#!/usr/bin/env bash
# Check status of Flower processes

echo "=== Flower Process Status ==="
echo ""
echo "SuperLink:"
ps aux | grep "flower-superlink" | grep -v grep || echo "  Not running"
echo ""
echo "SuperNodes:"
ps aux | grep "flower-supernode" | grep -v grep || echo "  Not running"
echo ""
echo "App Runner:"
ps aux | grep "flwr run" | grep -v grep || echo "  Not running"
echo ""
echo "=== Port Status ==="
netstat -tuln 2>/dev/null | grep -E "9092|9093|9094|9095" || echo "  No Flower ports in use"
