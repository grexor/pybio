# https://packaging.python.org/en/latest/tutorials/packaging-projects/
# https://stackoverflow.com/questions/50879668/python-setup-py-some-files-are-missing

rm dist/*
python -m build # build
pip install .
