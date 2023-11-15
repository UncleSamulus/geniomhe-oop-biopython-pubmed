
from Bio import Entrez
import matplotlib.pyplot as plt

def fetch_articles_by_keyword(keyword, start_date, end_date):
    Entrez.email = "your.email@example.com"  # Enter your email associated with NCBI
    Entrez.api_key = "your_api_key"  # Your API key from NCBI
    handle = Entrez.esearch(db="pubmed", term=keyword, mindate=start_date, maxdate=end_date, usehistory="y")
    record = Entrez.read(handle)
    handle.close()
    return record

def count_articles_by_time_interval(articles, interval):
    # Process article publication dates and count within the intervals
    # This part would involve organizing articles based on their publication dates

    # Return a dictionary or list with counts for each interval
    intervals = defaultdict(int)
    for article_id in articles['IdList']:
        summary = Entrez.esummary(db="pubmed", id=article_id)
        data = Entrez.read(summary)
        summary.close()
        pub_date = data['PubDate']
        year = pub_date.split()[0]  # Extracting the year from the publication date
        if interval == 'yearly':
            intervals[year] += 1
        # Add more intervals like monthly, quarterly, etc. if needed
    return intervals
def plot_trends(time_intervals, article_counts):
    # Use Matplotlib to create a plot of the trends
    # This might be a bar chart, line plot, etc. based on your chosen time intervals
    plt.bar(article_counts.keys(), article_counts.values())
    plt.xlabel('Year')
    plt.ylabel('Number of Articles')
    plt.title('Publication Trends Over Years')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()
# Usage
keyword = "your_keyword"
start_date = "20210101"  # Example start date in YYYYMMDD format
end_date = "20211231"    # Example end date in YYYYMMDD format

articles = fetch_articles_by_keyword(keyword, start_date, end_date)

# Choose a time interval (e.g., yearly, monthly) for counting articles
time_intervals = []  # Your time intervals
article_counts = count_articles_by_time_interval(articles, time_intervals)

# Plot the trends
plot_trends(time_intervals, article_counts)
