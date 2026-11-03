#!/usr/bin/env bash
echo "Installing required LaTeX packages..."
sudo tlmgr update --self
sudo tlmgr install $(cat tex_packages.txt)
echo "Done! You can now compile report.tex."