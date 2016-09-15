ENV_BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export PATH=$ENV_BASE_DIR/bin:$PATH
export PYTHONPATH=$ENV_BASE_DIR:$PYTHONPATH

# Avoid profiling in production situations
export PYTHONOPTIMIZE=2

unset ENV_BASE_DIR
