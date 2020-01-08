FROM python:3

ADD requirements.txt /requirements.txt
RUN set -x \
  && pip3 install -r /requirements.txt \
  && rm /requirements.txt

ADD app.py /app/app.py

ENV CREDENTIALS='{"user":"pass"}'
ENV HOST="0.0.0.0"
ENV DEBUG="false"

ENV PORT="80"
EXPOSE 80

ENV DATA="/data"
VOLUME [ "/data" ]

CMD [ "python3", "/app/app.py" ]
