#!/usr/bin/env bash
# Stop all Flower processes on this node

echo "Stopping Flower processes..."

# Stop SuperLink
pkill -f "flower-superlink" && echo "✓ Stopped SuperLink" || echo "No SuperLink process found"

# Stop SuperNodes
pkill -f "flower-supernode" && echo "✓ Stopped SuperNode" || echo "No SuperNode process found"

# Stop app runner
pkill -f "flwr run" && echo "✓ Stopped app runner" || echo "No app runner process found"

echo "Done. Verify with: ps aux | grep -E 'flower|flwr'"
