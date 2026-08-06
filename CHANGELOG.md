## [0.7.0] - 2026-08-06

### 🚀 Features

- Add ChromTable files parsing

### 🐛 Bug Fixes

- ChromTable field order
- Strands may be encoded as 1/-1
- Collision between strand indicator and landscape delimiter
- GFF attributes parsing may contain '='
- Correctly wrap IO errors
- Ignore empty GFF attributes

### 🚜 Refactor

- Avoid unwraps
- Remove useless variant

### ⚡ Performance

- Run the whole DB insertion in a single transaction

### ⚙️ Miscellaneous Tasks

- Update rusqlite
- Release syntesuite version 0.5.0
- Release syntesuite version 0.6.0
- Release syntesuite version 0.6.1
- Release syntesuite version 0.6.2
- Convert to devenv
- Update dependencies
- Clippy
## [0.4.0] - 2023-08-05

### 🚀 Features

- Add a way to access present species

### ⚙️ Miscellaneous Tasks

- Clippy
- Release syntesuite version 0.4.0
## [0.3.0] - 2023-07-16

### 🚀 Features

- Add strand information to the landscape
- Add strand, drop pretty landscape names

### ⚙️ Miscellaneous Tasks

- Release syntesuite version 0.3.0
## [0.2.4] - 2023-06-14

### ⚙️ Miscellaneous Tasks

- Downgrade dependencies for Guix
- Release syntesuite version 0.2.4
## [0.2.3] - 2023-06-14

### ⚙️ Miscellaneous Tasks

- Update dependencies
- Release syntesuite version 0.2.3
## [0.2.2] - 2023-05-29

### 🐛 Bug Fixes

- Obsolete column name

### ⚙️ Miscellaneous Tasks

- Add git-cliff as a dependency
- Release syntesuite version 0.2.2
## [0.2.1] - 2023-02-27

### 🐛 Bug Fixes

- Wrong column name

### ⚙️ Miscellaneous Tasks

- Release syntesuite version 0.2.1
## [0.2.0] - 2023-02-27

### 🚀 Features

- Add BED parsing
- Add more infos to genes

### 🚜 Refactor

- Share data structures

### ⚙️ Miscellaneous Tasks

- Add missing fields
- Release syntesuite version 0.2.0
## [0.1.0] - 2023-01-23

### 🚀 Features

- Add dbmaker
- Add GeneBook
- ID column is customizeable

### ⚙️ Miscellaneous Tasks

- Setup everything
- Downgrade to Guix
- Release syntesuite version 0.1.0
