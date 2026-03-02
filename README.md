# User experience

## What is sgevalviz?
An improved version of sgeval with visualization and statistical analysis tools for GTFS bioinformatics data.
You can find the official PyPi [here](https://pypi.org/project/sgevalviz/), and the Test PyPi [here](https://test.pypi.org/project/sgevalviztest/).

## How to execute
To execute the program you must run: `sgevalviz <path/to/the/folder/you/want/to/store/the/data> <path/to/candidate/file.gtf> <path/to/baseline/file.gtf>`. If you are using `sgevalviztest` you run the command above the same but substitute `sgevalviz` by `sgevalviztest`.
Extra parameters can be add to the command, for example `--candidate-config=augustus` and it will use the config defined on `src/sgevalviz/configs/augustus.json` for the `<path/to/candidate/file.gtf>`, you can substiture `candidate` by `baseline` to affect the baseline data instead. You can also pass the full path of your custom config with `--custom-candidate-config=~/configs/augustus.json`, for example, and substitute by `baseline` if needed.


# Developer experience
To develop you must access the Github shown in the PyPi description and be able to commit and push. Only Linux and Mac are tested to work with the scripts bellow.

## How to test
If you are developing locally, you can test the code by runining `./development.sh <path/to/the/folder/you/want/to/store/the/data> <path/to/candidate/file.gtf> <path/to/baseline/file.gtf>`. It's the same as with the official or test one just substituting the command by `./development.sh`.

## How to publish to PyPi or Test PyPi
After you already tested locally and are done commiting and pushing, you'll need to change versions and tag it so the workflow can publish it to PyPi or Test PyPi. If the last version on `pyproject.toml` is official (e.g. the version is "A.B.C" and not "A.B.CrcD" where A, B, C and D are numbers), you will need to increase the major, minor or patch version. To do this, you'll run `./bump_version.sh <upgrade_name>` where `upgrade_name` is either `major`, `minor` or `patch`. After this, you'll increase the version and will append `rc1` to it, this will make the code being pushed to the Github and the workflow will make the code available at Test PyPi - you can find the links above. You should then test the code by updating `sgevalviztest` that you can find on Test PyPi, and if everything works fine you can tag it as an official version, otherwise, you'll make the needed changes and bump the `D` in the `rcD` value. To do this, you'll execute `./bump_version.sh <upgrade_name>` but now `<upgrade_name>` must be either `rc` or `official`. `rc` will increase the `D` in the `rcD` value and update the Test PyPi version and `official` will remove the `rcD` ending, making it official and updating the official PyPi version.
