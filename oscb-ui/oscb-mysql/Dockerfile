FROM mysql:8.0
ENV MYSQL_ROOT_PASSWORD Shirley@2410
COPY oscb_schema_latest.sql /docker-entrypoint-initdb.d/
# CMD ["mysqld"]