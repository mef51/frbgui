Build a preview with `make html && open _build/html/index.html`

Use `livereloadx -s _build/html` in one terminal to start a server and `make html` in the other so you don't have a new tab every change.

Do `make index` to regenerate the index and make the html.

Read the docs does this command:
```
python -m sphinx -T -W --keep-going -b html -d _build/doctrees -D language=en . docs/_build/html
```
