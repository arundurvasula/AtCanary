#AtCanary

Using dadi to create a demographic model of the introduction of Arabidopsis to the Canary islands.

Organization:

data: contains the data!
results: contains results!
scripts:
  - `Arabidopsis.ipynb`: IPython notebook containing code to test out models and explain what we're doing (may contain pretty pictures)
  - `demographic_models.py`: python script that contains the model functions we are working with. Should contain all of our models.
  - `simple_unfolded.py`: a script to run the simple canary island model.


# Simple Canary Island Model
This model is just a 3 population split with no migration. Instantaneous population size changes.

```
    |
     ----
    |    |
    |     ----
    |    |    |
    S    M    C
```

# Simple Canary Island Model with migration
This model is just a 3 population split with migration. Instantaneous population size changes.

```
    |
     ----
    |    |
    |     ----
    |====|===>|
    S    M    C
```
