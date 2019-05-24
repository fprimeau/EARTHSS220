This ReadMe is intended for the mainainers of the repository.

This repository is structured in the following way:

```
EARTHSS220
├── ReadMe.md
└── src
    ├── ReadMe.md
    ├── generate.jl
    ├── lectures
    │   ├── lecture_0
    │   │   └── prerequisites.md
    │   └── lecture_1
    │       └── ideal_mean_age.jl
    └── generated
        ├── lecture_0
        │   └── prerequisites.md
        └── lecture_1
            └── ideal_mean_age.ipynb
```

You should only modify the source files in `src/lectures/`.

To generate the notebooks from the Literate-ready `.jl` files in `src/lectures/`, from the root of the repository (in `EARTHSS220/`), run the `generate.jl` script via

```bash
julia src/generate.jl
```

Adding the links to the notebooks in the main `ReadMe.md` is done manually.