import re
import numpy as np
import pandas as pd


def remove_accents(string: str) -> str:
    
    string = re.sub(r'[ÁáÀà]', 'A', string)
    string = re.sub(r'[Éé]', 'E', string)
    string = re.sub(r'[Íí]', 'I', string)
    string = re.sub(r'[Óó]', 'O', string)
    string = re.sub(r'[Úú]', 'U', string)

    return string

def swap_words(string: str, split_pattern: str = '. ') -> str:
    
    words = string.split(split_pattern)
    if len(words) == 2:
        return ' '.join([word.strip() for word in words[::-1]])
    else:
        return string

def arabic_to_roman(match):
    arabic = int(match.group(0))
    roman_numerals = {
        1: 'I', 4: 'IV', 5: 'V', 9: 'IX', 10: 'X', 40: 'XL',
        50: 'L', 90: 'XC', 100: 'C', 400: 'CD', 500: 'D', 900: 'CM', 1000: 'M'
    }
    result = ''
    for value, numeral in sorted(roman_numerals.items(), key=lambda x: -x[0]):
        while arabic >= value:
            result += numeral
            arabic -= value
    return result

def correct_names(df: pd.DataFrame, col_pattern: str = 'nombre', split_pattern: str = ', ') -> pd.DataFrame:
    
    col_names = [col for col in df if col_pattern in col.lower()]
    for col in col_names:
        # replace missing values
        df[col] = df[col].replace(np.nan, '')
        # remove accents
        df[col] = df[col].apply(remove_accents)
        # swap articles
        df[col] = df[col].apply(swap_words, split_pattern=split_pattern)
        # convert arabic numbers to roman numbers
        df[col] = df[col].apply(lambda x: re.sub(r'\b\d+\b', arabic_to_roman, x))
        
    return df