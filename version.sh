printf 'module recontig._version;\nenum VERSION="%s";\n' $(git describe --tags) > source/recontig/version.d
