## timelineData

- Getting a lot of data still spends quite a while in GC (20% over 170,000 items)
  - It is still reasonably fast, less than 3 seconds
- It would be useful to not require Float64, removing this creates some type instability
  - I think this needs a type on DataMatrix and MarketData
- Once dropmissing is called, a lot of the functions stop working correctly
  - dropmissing might need to return a different type, perhaps AxisArrays?
  - Storage and access might just need to be two different types...
  - Looking closer at MatrixTable which is default in Tables.jl might make the most sense
- I also wonder about building this on top of AxisArrays
  - My only concern is joining and missing data
    - Doing a test hcat, it returns a simple Matrix, which might not fit my needs since I would need to build a table interface on top of that

## fastRegression

- If there are too few observations, I think it is best to return an object
  - The object would have number of observations and all else undefined
  - This prevents needing functions that accept a missing type
    - Though these functions would need a check on undefined and return missing instead
