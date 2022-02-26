## timelineData

- I think what would help is to make the inital calls more lazy
  - StatsModels just needs access to a vector, so no need to build a matrix multiple times
- Store the data in Dicts that access a single column
- When the data is accessed, do not copy the data, just reference the parent and the data
  - This will need to be a mutable struct to allow for dropmissing (which is just a Bool)
  - also need some column selection

## fastRegression

- If there are too few observations, I think it is best to return an object
  - The object would have number of observations and all else undefined
  - This prevents needing functions that accept a missing type
    - Though these functions would need a check on undefined and return missing instead
