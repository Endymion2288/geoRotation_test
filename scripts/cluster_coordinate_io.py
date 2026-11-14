"""Read the frozen cluster audit, including its recorded Gaudi label encoding."""
import gzip
import json


def read_coordinates(path):
    with gzip.open(path, 'rt') as file:
        for line in file:
            # Gaudi::Property<double>::operator<< included its name in v1.
            # Remove only that literal label; the actual numerical token is kept.
            line = line.replace('"variance_scale": \'MeasurementVarianceScale\':',
                                '"variance_scale":')
            yield json.loads(line)
