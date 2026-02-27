#!/bin/bash
set -e

INSTALL_DIR="$HOME/.local/bin"
ENV_NAME="gatk3"

for arg in "$@"; do
	case "$arg" in
		--install-dir=*) INSTALL_DIR="${arg#*=}" ;;
		--env-name=*)    ENV_NAME="${arg#*=}" ;;
		*) echo "Unknown argument: $arg"; exit 1 ;;
	esac
done

mkdir -p "$INSTALL_DIR"

echo "Creating conda environment '$ENV_NAME' from pseudoDB_env.yaml..."
conda env create -n "$ENV_NAME" -f pseudoDB_env.yaml

echo "Installing GATK 3.8 into environment '$ENV_NAME'..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
tar -jxvf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
mv GenomeAnalysisTK-3.8-1-0-gf15c1c3ef gatk3.8
gatk-register gatk3.8/GenomeAnalysisTK.jar
rm -rf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 gatk3.8
conda deactivate

echo "Installing pseudoDB to $INSTALL_DIR..."
SOURCE_DIR="$(cd "$(dirname "$0")" && pwd)"
cat > "$INSTALL_DIR/pseudoDB" <<EOF
#!/bin/bash
python3 "$SOURCE_DIR/script/pseudoDB.py" "\$@"
EOF
chmod +x "$INSTALL_DIR/pseudoDB"

echo "Installation complete."
echo "Please activate the environment before use: conda activate $ENV_NAME"

if [[ ":$PATH:" != *":$INSTALL_DIR:"* ]]; then
	echo "Warning: $INSTALL_DIR is not on your PATH."
	echo "Add this to your ~/.bashrc or ~/.zshrc:"
	echo "  export PATH=\"$INSTALL_DIR:\$PATH\""
fi
