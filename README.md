# mcc: Modal Clustering for Categorical Data

**mcc** implements a **novel clustering approach for categorical data**. Traditional clustering methods often struggle with categorical variables, especially when the number of clusters is unknown. Our approach introduces a **new notion of clusters** based on two complementary principles:

1. **High frequency**
2. **Variable association**

This dual concept aligns with the **modal clustering framework** used in continuous data, allowing us to adapt existing operational tools to the categorical setting.  

### Key Features

- Automatically identifies the **number of clusters**—no need to pre-specify it.  
- Works directly with **categorical data** without requiring numeric encoding.  
- Leverages **mode-shift and level-set approaches** to detect clusters robustly.
- **Fast and efficient** even with **large datasets**

The method is suitable for datasets with a stong dependency structure, a high association between the variables or an unknown cluster structure.  

---

## Installation

You can install the package directly from GitHub using the `remotes` package:

```r
# Install devtools if you don't have it
install.packages("devtools")

# Install mcc from GitHub
library(devtools)
install_github("NoeCors/mcc")
```

## Usage
```r
library(mcc)
data(titanic)
head(titanic)

### Apply modal clustering
clusters <- mcc(titanic)
```

## References

Corsini, N. & Menardi, G. (2025). *Modal clustering for categorical data*.
