import requests
import json
from typing import Tuple
import pandas as pd
import datetime as dt
import pytz
import os

def save_metadata(
                data, 
                path, 
                name
                ):
    path = '/'.join([path, f'{name}.json'])
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=4, ensure_ascii=False)

def yahoo_scraper(
                ticker_, 
                interval_ = '5m', 
                range_ = '60d'
                ) -> Tuple[ pd.DataFrame, dict, str, str]:
    """Função para realizar web scraping do Yahoo Finance de um ticker (ativo) específico.

    Args:
        ticker_ (_type_): identificador do ativo na bolsa
        interval_ (str, optional): frequência dos dados. padrão de 5-5 minutos (5m).
        range_ (str, optional): período dos dados. padrão de últimos 60 dias (60d).

    Returns:
        Tuple[str, pd.DataFrame, str]: tupla de 4 valores com o dataframe gerado, metadados do ativo, 
        o título do ativo sem caracteres especiais e a hora de extração, respectivamente.
    """

    # Especificações de página
    url = f'https://query1.finance.yahoo.com/v8/finance/chart/{ticker_}?interval={interval_}&range={range_}'
    headers = {
        'User-Agent': 'Mozilla/5.0'
    }

    print(f"\t Extraindo a série {ticker_} com frequência {interval_} para o período {range_}...")

    # Dados da página
    extraction_time = dt.datetime.now().strftime("%Y%m%d%H")

    try:
        response = requests.get(url, headers=headers)
        data = response.json()
        meta = data['chart']['result'][0]['meta']
        timestamp = data['chart']['result'][0]['timestamp']
        indicators = data['chart']['result'][0]['indicators']['quote'][0]
    except Exception as e:
        print(f"\t\t Problema de extração: {e}")
        return


    # Estruturando os dados obtidos
    df = pd.DataFrame(indicators)

    df['datetime_utc'] = pd.to_datetime(timestamp, unit='s', utc=True)
    sao_paulo = pytz.timezone("America/Sao_Paulo") # ajuste para janela de negociacao do Brasil
    df['datetime'] = df['datetime_utc'].dt.tz_convert(sao_paulo)
    df['datetime'] = df['datetime'].dt.tz_localize(None)
    df.set_index('datetime', inplace=True)
    
    # Identificadores da série
    ticker_formatted = ticker_.replace('^', '').replace('.','').replace('-','_')
    df['ticker'] = ticker_formatted
    df['frequency'] = interval_
    df['data_extracao'] = extraction_time

    # Seleciona colunas
    columns = ['ticker', 'frequency', 'open', 'high', 'low', 'close', 'data_extracao']
    df = df.loc[:,columns]

    return df, meta, ticker_formatted, extraction_time 

def list_tickers() -> dict:
    # Variáveis de interesse
    d = {
            '^BVSP',           # Ibovespa
            'BOVA11.SA',       # iShares Ibovespa ETF
            '^GSPC',           # S&P 500 (for comparison)
            'IFIX.SA',
            'SMLL.SA',
            'IDIV.SA',

            # Mining & Oil
            'VALE3.SA',        # Vale S.A. (iron ore, nickel)
            'PETR3.SA',        # Petrobras ON (oil)
            'PETR4.SA',        # Petrobras PN (oil)

            # Pulp & Paper
            'SUZB3.SA',        # Suzano (pulp)
            'KLBN11.SA',       # Klabin (pulp and paper)

            # Steel & Mining
            'CSNA3.SA',        # CSN (iron, steel)
            'GGBR4.SA',        # Gerdau (steel)

            # Agribusiness (meat, grains, bioenergy)
            'JBSS3.SA',        # JBS (meat, grain input)
            'BRFS3.SA',        # BRF (poultry/pork, grain input)
            'BEEF3.SA',        # Minerva (beef, grain input)
            'RAIZ4.SA',        # Raízen (sugarcane, ethanol)

            # Grain producers
            'SLC3.SA',         # SLC Agrícola (soy, corn, cotton)
            'AGRO3.SA',         # BrasilAgro (grains, land)

            # Crypto
            "BTC",  # Bitcoin
            "ETH"  # Ethereum
        }
    return d

###############################################################################
########### SCRAPER ###########

# Variáveis de interesse
tickers = list_tickers()

# Frequência do dado + intervalo disponível na origem
intervals = {
            '5m' : '60d'
            }

print('Iniciando extrações do Yahoo Finance.')

# Realiza download das séries especificadas
path = 'data/raw/'
for ticker in tickers:
    for interval in intervals.keys():
        try:
            df, meta, ticker_formatted, extraction_time = yahoo_scraper(ticker, interval, intervals[interval])
        except Exception as e:
            continue
        save_metadata(meta, path, f"metadata_{ticker_formatted}")
        title = f'{path}/{ticker_formatted}_{interval}_{extraction_time}'
        df.to_parquet(f'{title}.parquet')

print('Extrações concluídas.')

#################################################################################
########### TRATAMENTO ###########
# Concatena dados
in_path = 'data/raw/'
out_path = 'data/processed/'

# df = pd.DataFrame()
df = pd.read_parquet('data/raw/history.parquet').reset_index()
tickers = list_tickers()
for ticker in tickers:
    ticker_str = ticker.replace('^', '').replace('.','').replace('-','_')
    ticker_hist = pd.DataFrame()
    ticker_files = [x for x in os.listdir(in_path) if ticker_str in x and '.parquet' in x]
    for ticker_file in ticker_files:
        filepath = '\\'.join([in_path, ticker_file])
        df_ticker = pd.read_parquet(filepath)
        df_ticker.reset_index(inplace=True)
        ticker_hist = pd.concat([ticker_hist, df_ticker])
    ticker_hist.drop_duplicates(inplace=True)
    ticker_hist.to_excel(f'{out_path}{ticker_str}.xlsx', index=False)
    df = pd.concat([df, ticker_hist], axis=0, ignore_index=True)

# df.to_excel(f'{out_path}database.xlsx', index=False)
df.to_parquet(f'{out_path}database.parquet', index=False)

print('Dados agregados.')
