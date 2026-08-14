## [0.8.0] - 2026-08-14

### 🚜 Refactor

- Strongly type family IDs
## [0.7.0] - 2026-08-06

### 🐛 Bug Fixes

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

- Convert to devenv
- Update dependencies
- Clippy
- Release syntesuite version 0.7.0
## [0.6.2] - 2024-08-05

### 🐛 Bug Fixes

- Strands may be encoded as 1/-1

### ⚙️ Miscellaneous Tasks

- Release syntesuite version 0.6.2
## [0.6.1] - 2024-07-30

### 🐛 Bug Fixes

- ChromTable field order

### ⚙️ Miscellaneous Tasks

- Release syntesuite version 0.6.1
## [0.6.0] - 2024-07-30

### 🚀 Features

- Add ChromTable files parsing

### ⚙️ Miscellaneous Tasks

- Release syntesuite version 0.6.0
## [0.5.0] - 2024-07-21

### ⚙️ Miscellaneous Tasks

- Update rusqlite
- Release syntesuite version 0.5.0
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
