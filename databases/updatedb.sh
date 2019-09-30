#!/bin/bash
mongo 127.0.0.1:27017/rylanTableV1 rylanClassifierDb.js
mongo 127.0.0.1:27017/rylanFiltersV1 filtersTable.js
mongo 127.0.0.1:27017/rylanTableV1 rylanClassifierInDelDb.js
mongo 127.0.0.1:27017/rylanFiltersV1 filtersTableInDel.js
